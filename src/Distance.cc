#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>

#include "Distance.h"
#include "Seq.h"

using namespace std;

Distance::Distance(int kmer, double threshold, int step_size) {
    this->k = kmer;
    this->thrs = threshold;
    this->step = step_size;
}

bool Distance::compare(const Seq& s, const Seq& t) {
    return distance(s, t) >= thrs;
}

double Distance::distance(const Seq& s, const Seq& t) {
    size_t slen = s.length,
           tlen = t.length;

    size_t long_len, short_len;
    if (slen >= tlen) {
        long_len = slen;
        short_len = tlen;
    } else {
        long_len = tlen;
        short_len = slen;
    }

    uint8_t *longer  = slen >= tlen ? s.data : t.data;
    uint8_t *shorter = slen >= tlen ? t.data : s.data;

    // allocate array of length equal to the number of different kmers
    static const int kmer_count = pow(4, k);
    int *kmers = new int[kmer_count](); // zero initialized due to ()

    static const uint32_t k2 = 2 * k;
    static const uint32_t mask = kmer_count - 1;    // 0b001111 (2*k 1's)
    int cur_dist = 0;

    // count kmers in the shorter and the longer string, respectively
    for (size_t i = 0; i <= short_len - k; ++i) {
        uint32_t kmer_l = 0; // binary repr. of kmer in longer sequence
        uint32_t kmer_s = 0;

        //memcpy(&kmer_l, (longer + (i/4)), 4);
        kmer_l = stream2int(longer + (i/4));
        kmer_l >>= (32 - 2*(i % 4) - k2);
        kmer_l &= mask;

        kmer_s = stream2int(shorter + (i/4));
        kmer_s >>= (32 - 2*(i % 4) - k2);
        kmer_s &= mask;

        kmers[kmer_l]++ < 0 ? --cur_dist : ++cur_dist;
        kmers[kmer_s]-- > 0 ? --cur_dist : ++cur_dist;

    }

/*    cout << endl << cur_dist << " : ";
    cur_dist = 0;
    // Manhattan distance between the two strings
    for (int i = 0; i < kmer_count; ++i)
        cur_dist += abs(kmers[i]);*/

    int min_dist = cur_dist;    // the least distance window so far
    int win_size = short_len;
    int windows = long_len - short_len;
    int total = 2 * (short_len - k + 1);

    double jaccard_dist = 0;

    if (windows == 0) {
        delete[] kmers;
        return (double) (total - cur_dist) / (double) total;
    }

    /*
     * pre_gram:  kmer moving out of window
     * post_gram: kmer moving into window
     *
     * actgactgactg
     * actgactgactgactgactg
     * ^^^^     ^^^^
     * pre      post
     */
    for (int i = 0; i < windows; ++i) {
        uint32_t pre_gram = 0;
        uint32_t post_gram = 0;
        int post_i = i + win_size - k + 1;

        pre_gram = stream2int(longer + i/4);
        pre_gram >>= (32 - 2*(i % 4) - k2);
        pre_gram &= mask;

        post_gram = stream2int(longer + post_i/4);
        post_gram >>= (32 - 2*(post_i % 4) - k2);
        post_gram &= mask;

        if (pre_gram == post_gram)
            continue;   // same kmers, so no need to calculate new distance

        // if changed for the better, decrement cur_dist, otherwise increment
        kmers[pre_gram]  > 0 ? --cur_dist : ++cur_dist;
        kmers[post_gram] < 0 ? --cur_dist : ++cur_dist;

        // adjust kmer count from change
        --kmers[pre_gram];
        ++kmers[post_gram];

        min_dist = min(cur_dist, min_dist);

        jaccard_dist = (double) (total - min_dist) / (double) total;
        if (jaccard_dist >= thrs) {
            delete[] kmers;
            return jaccard_dist;
        }
    }
    delete[] kmers;
    return jaccard_dist;
}

/* Returns a sorted vector by decreasing order and returns the n most
   frequent kmers if they exist */
vector<int> Distance::compute_key(const Seq& s, int n) {
    const size_t slen = s.length;
    uint8_t *data = s.data;

    // allocate array of length equal to the number of different kmers
    const int kmer_count = pow(4, k);
    //int *kmers = new int[kmer_count](); // zero initialized due to ()
    vector<int> kmers;
    kmers.resize(kmer_count);

    static const uint32_t k2 = 2 * k;
    static const uint32_t mask = pow(2, k2) - 1; // 0b001111 (2*k 1's)

    // count kmers in the sequence
    for (size_t i = 0; i <= slen - k; ++i) {
        uint32_t kmer = 0; // binary repr. of kmer

        kmer = stream2int(data + (i/4));
        kmer >>= (32 - 2*(i % 4) - k2);
        kmer &= mask;

        ++kmers[kmer];
    }

    vector<int> ret;
    ret.resize(n);

    partial_sort_copy(kmers.begin(), kmers.end(),
            ret.begin(), ret.end(), greater<int>());


    //delete[] kmers;
    return ret;
}

double Distance::levenshtein(const string& s, const string& t) {
    int slen = s.length();
    int tlen = t.length();

    // Trivial cases
    if (slen == 0) {
        return tlen;
    }
    if (tlen == 0) {
        return slen;
    }
    int col[slen + 1];
    int pcol[slen + 1];
    for (int i = 0; i < slen; i++) {
        pcol[i] = i;
    }
    // Dynamic approach to calculate the distance between two strings
    for (int i = 0; i < tlen; i++) {
        col[0] = i+1;
        for (int j = 0; j < slen; j++) {
            int cost = ((s[j]) == (t[i])) ? 0 : 1;
            col[j+1] = min(col[j] + 1, min(pcol[j+1] + 1, pcol[j] + cost));
        }
        for (int j = 0; j < slen + 1; j++) {
            pcol[j] = col[j];
        }
    }
    if (tlen > slen)
        return (double)(tlen - col[slen]) / (double)tlen;
    else
        return (double)(slen - col[slen]) / (double)slen;
    //return col[slen];
}

double Distance::levenshtein_window(const string& s, const string& t) {
    int slen = s.length();
    int tlen = t.length();

    string shorter, longer;
    int short_len, long_len;
    double min_dist = 9999;
    double cur_dist;

    if (slen <= tlen) {
        shorter = s;
        short_len = slen;
        longer = t;
        long_len = tlen;
    } else {
        shorter = t;
        short_len = tlen;
        longer = s;
        long_len = slen;
    }

    int win_size = short_len;
    int windows = long_len - short_len;

    if (windows == 0) {
        return (double)(win_size - levenshtein(s, t)) / (double)win_size;
    }

    for (int i = 0; i < windows; i++) {
        cur_dist = levenshtein(shorter, longer.substr(i, win_size));
        min_dist = min(min_dist, cur_dist);
    }

    return (double)(win_size - min_dist) / (double)win_size;
}

/*void Distance::printDistMatrix(const char* filename, int count) {
    ifstream fs0(filename);
    ifstream fs1(filename);

    string fst, snd;
    double distances[count][count];

    for (int i = 0; i < count; i++)
        for (int j = 0; j < count; j++)
            distances[i][j] = -1; // initialize distance matrix entries to -1

    for (int i = 0; i < count; ++i) {
        io::read_sequence(fs0, fst);
        transform(fst.begin(), fst.end(), fst.begin(), ::tolower);

        for (int j = 0; j < count; ++j) {
            io::read_sequence(fs1, snd);
            transform(snd.begin(), snd.end(), snd.begin(), ::tolower);
            if (i == j) { // don't compare a sequence to itself
                distances[i][j] = 0;
                continue;
            }
            if (distances[i][j] != -1) {
                continue;
            }

            double newdist = levenshtein(fst, snd);
            distances[i][j] = newdist;
            distances[j][i] = newdist;
        }
        fs1.seekg(0, ios::beg); // rewind fs1 to start
    }

    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < count; ++j) {
            cout << " " << distances[i][j]; // TODO: reset width?
        }
        cout << '\n';
    }

    fs0.close();
    fs1.close();
}

void Distance::jac_printDistMatrix(const char* filename, int count) {
    ifstream fs0(filename);

    double distances[count][count];
    vector<Seq> seqs;

    io::read_seqs(fs0, seqs, count);

    for (int i = 0; i < count; i++)
        for (int j = 0; j < count; j++)
            distances[i][j] = -1; // initialize distance matrix entries to -1

    for (int i = 0; i < count; ++i) {

        for (int j = 0; j < count; ++j) {
            if (i == j) { // don't compare a sequence to itself
                distances[i][j] = 0;
                continue;
            }
            if (distances[i][j] != -1) {
                continue;
            }

            double newdist = distance(seqs[i], seqs[j]);
            distances[i][j] = newdist;
            distances[j][i] = newdist;
        }
    }

    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < count; ++j) {
            cout << " " << distances[i][j]; // TODO: reset width?
        }
        cout << '\n';
    }

    fs0.close();
}*/

/*set<string> Distance::kmers(const Seq& s) {
    set<string> kmers;
    for (unsigned int i = 0; i <= s.data.length()-k; i++) {
        kmers.insert(s.data.substr(i,k));
    }
    return kmers;
}*/
