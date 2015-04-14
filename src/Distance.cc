#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>
#include <bitset>

#include "Distance.h"
#include "IO.h"

using namespace std;

Distance::Distance(int kmer, double threshold, int step_size) {
    this->k = kmer;
    this->thrs = threshold;
    this->step = step_size;
}

inline uint8_t nth_left_2bits(uint8_t b, const int& n) {
    int shift = 6 - 2 * (n % 4);
    return b & (3 << shift) >> shift;
}

bool Distance::compare(const Seq& s, const Seq& t) {
    size_t slen = s.length(),
           tlen = t.length();

    size_t long_len, short_len;
    if (slen >= tlen) {
        long_len = slen;
        short_len = tlen;
    } else {
        long_len = tlen;
        short_len = slen;
    }

    uint8_t *longer  = slen >= tlen ? s.data() : t.data();
    uint8_t *shorter = slen >= tlen ? t.data() : s.data();

    // allocate array of length equal to the number of different kmers
    const int kmer_count = pow(4, k);
    int *kmers = new int[kmer_count](); // zero initialized due to ()
    //unordered_map<uint32_t, int> kmers;
    /*vector<int> kmers = new vector<int>;
    kmers.resize(pow(4,k));*/

    // count kmers in the shorter and the longer string, respectively
    for (size_t i = 0; i <= short_len - k; ++i) {
        uint32_t kmer_l = 0; // binary repr. of kmer in longer sequence
        uint32_t kmer_s = 0;

        //kmer_l |= (*longer  & (3 << shift) >> shift;
        //kmer_s |= (*shorter & (3 << shift) >> shift;
        for (int j = 0; j < k; ++j) {
            int shift = 6 - 2 * ((i + j) % 4);
            kmer_l |= (*(longer  + (i+j)/4) & (3 << shift)) >> shift;
            kmer_s |= (*(shorter + (i+j)/4) & (3 << shift)) >> shift;
            if (j == k-1)
                break;
            kmer_l <<= 2;
            kmer_s <<= 2;
        }

        ++kmers[kmer_l];
        --kmers[kmer_s];
    }

    // Manhattan distance between the two strings
    /*int cur_dist = 0;
    for (auto it = kmers.cbegin(); it != kmers.cend(); ++it)
        cur_dist += abs(it->second);
        //cur_dist += abs(*it);*/
    int cur_dist = 0;
    for (int i = 0; i < kmer_count; ++i)
        cur_dist += abs(kmers[i]);

    int min_dist = cur_dist;    // the least distance window so far
    int win_size = short_len;
    int windows = long_len - short_len;
    int total = 2 * (short_len - k + 1);

    if (windows == 0) {
        delete[] kmers;
        return ((double) (total - cur_dist) / (double) total) >= thrs;
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
        for (int j = 0; j < k; ++j) {
            int shift = 6 - 2 * ((i + j) % 4);
            int post_i = i + win_size - k + 1;
            int post_shift = 6 - 2 * ((post_i + j) % 4);
            pre_gram |= (*(longer + (i+j)/4) & (3 << shift)) >> shift;
            post_gram |= (*(longer + (post_i+j)/4) & (3 << post_shift)) >> post_shift;
            if (j == k-1)
                break;
            pre_gram <<= 2;            
            post_gram <<= 2;
        }

        if (pre_gram == post_gram)
            continue;   // same kmers, so no need to calculate new distance

        // if changed for the better, decrement cur_dist, otherwise increment
        kmers[pre_gram]  < 0 ? --cur_dist : ++cur_dist;
        kmers[post_gram] > 0 ? --cur_dist : ++cur_dist;
    
        // adjust kmer count from change
        ++kmers[pre_gram];
        --kmers[post_gram];

        min_dist = min(cur_dist, min_dist);
        cout << cur_dist << endl;

        if (((double) (total - min_dist) / (double) total) >= thrs) {
            delete[] kmers;
            return true;
        }
    }

    delete[] kmers;
    return false;
}

/**
 * Given a string, return an int representation of the string,
 * e.g. gram_pos("acgt") == 0b11100100. 'a' == 0b00, 'c' == 0b01 etc.
 */
inline unsigned int gram_pos(const string& s) {
    int slen = s.length();
    unsigned int res = 0;

    for (int i = 0; i < slen; ++i) {
        switch (s[i]) {
            case 'c': case 'C':
                res |= 1 << (i*2);  // 0b01
                break;
            case 'g': case 'G':
                res |= 2 << (i*2);  // 0b10
                break;
            case 't': case 'T': case 'u': case 'U':
                res |= 3 << (i*2);  // 0b11
                break;
            default:
                break;
        }
    }
    return res;
}

/**
 * Given two strings, the int k in k-mer (word length), calculated the
 * d2-distance between the strings based on k-mers: count the occurences of
 * each poosible k-mer, calculate the Manhattan distance between the two k-mer
 * occurence vectors. Return true if within threshold in some window (of size
 * shortest string), otherwise return false.
 */
bool Distance::compare(const vector<bitset<2>>& s, const vector<bitset<2>>& t) {
    int slen = s.size(),
        tlen = t.size();

    int long_len, short_len;
    if (slen >= tlen) {
        long_len = slen;
        short_len = tlen;
    } else {
        long_len = tlen;
        short_len = slen;
    }

    const vector<bitset<2>>& longer = slen >= tlen ? s : t;

    // allocate array of length equal to the number of different kmers
    const int kmer_count = pow(4, k);
    int *kmers = new int[kmer_count](); // zero initialized due to ()

    // count kmers in the shorter and the longer string, respectively
    auto it_s = s.cbegin(), it_t = t.cbegin();
    for (int i = 0; i <= short_len-k; ++i) {
        unsigned long pos_s = (*it_s).to_ulong();
        unsigned long pos_t = (*it_t).to_ulong();
        for (int i = 0; i < k-1; ++i) {
            pos_s <<= 2;
            pos_s |= (*(it_s+i)).to_ulong();
            pos_t <<= 2;
            pos_t |= (*(it_t+i)).to_ulong();
        }
        ++kmers[pos_s];
        --kmers[pos_t];
        ++it_s;
        ++it_t;
    }

    // count kmers in the shorter and the longer string, respectively
    /*for (int i = 0; i <= short_len-k; ++i) {
        ++kmers[gram_pos(shorter.substr(i,k))];
        --kmers[gram_pos(longer.substr(i,k))];
    }*/

    // Manhattan distance between the two strings
    int cur_dist = 0;
    for (int i = 0; i < kmer_count; ++i)
        cur_dist += abs(kmers[i]);

    int min_dist = cur_dist;    // the least distance window so far
    int win_size = short_len;
    int windows = long_len - short_len;
    int total = 2 * (short_len - k + 1);

    if (windows == 0) {
        delete[] kmers;
        return cur_dist >= thrs;
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
    unsigned int post_pos;
    unsigned long pre_gram, post_gram;
    for (int i = 0; i < windows; ++i) {
        post_pos = i + win_size - k + 1;
        pre_gram = longer[i].to_ulong();
        post_gram = longer[post_pos].to_ulong();
        for (int j = 0; j < k-1; ++j) {
            pre_gram <<= 2;
            pre_gram |= longer[i+j].to_ulong();
            post_gram <<= 2;
            post_gram |= longer[post_pos+j].to_ulong();
        }

        /*pre_gram  = gram_pos(longer.substr(i, k));
        post_gram = gram_pos(longer.substr(i + win_size - k + 1, k));*/

        if (pre_gram == post_gram)
            continue;   // same kmers, so no need to calculate new distance

        // if changed for the better, decrement cur_dist, otherwise increment
        kmers[pre_gram]  < 0 ? --cur_dist : ++cur_dist;
        kmers[post_gram] > 0 ? --cur_dist : ++cur_dist;

        // adjust kmer count from change
        ++kmers[pre_gram];
        --kmers[post_gram];

        min_dist = min(cur_dist, min_dist);

        if (((double) (total - min_dist) / (double) total) >= thrs) {
            delete[] kmers;
            return true;
        }
    }

    delete[] kmers;
    return false;
}


/**
 * Given two strings, the int k in k-mer (word length), calculated the
 * d2-distance between the strings based on k-mers: count the occurences of
 * each poosible k-mer, calculate the Manhattan distance between the two k-mer
 * occurence vectors. Return true if within threshold in some window (of size
 * shortest string), otherwise return false.
 */
bool Distance::compare(const string& s, const string& t) {
    int slen = s.length(),
        tlen = t.length();
    string shorter, longer;
    int short_len, long_len;

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

    // allocate array of length equal to the number of different kmers
    int kmer_count = pow(4, k);
    int *kmers = new int[kmer_count](); // zero initialized due to ()

    // count kmers in the shorter and the longer string, respectively
    for (int i = 0; i <= short_len-k; ++i) {
        ++kmers[gram_pos(shorter.substr(i,k))];
        --kmers[gram_pos(longer.substr(i,k))];
    }

    // Manhattan distance between the two strings
    int cur_dist = 0;
    for (int i = 0; i < kmer_count; ++i)
        cur_dist += abs(kmers[i]);

    int min_dist = cur_dist;    // the least distance window so far
    int win_size = short_len;
    int windows = long_len - short_len;
    int total = 2 * (short_len - k + 1);

    if (windows == 0) {
        delete[] kmers;
        return cur_dist >= thrs;
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
    unsigned int pre_gram, post_gram;
    for (int i = 0; i < windows; ++i) {
        pre_gram  = gram_pos(longer.substr(i, k));
        post_gram = gram_pos(longer.substr(i + win_size - k + 1, k));

        if (pre_gram == post_gram)
            continue;   // same kmers, so no need to calculate new distance

        // if changed for the better, decrement cur_dist, otherwise increment
        kmers[pre_gram]  < 0 ? --cur_dist : ++cur_dist;
        kmers[post_gram] > 0 ? --cur_dist : ++cur_dist;
    
        // adjust kmer count from change
        ++kmers[pre_gram];
        --kmers[post_gram];

        min_dist = min(cur_dist, min_dist);

        if (((double) (total - min_dist) / (double) total) >= thrs) {
            delete[] kmers;
            return true;
        }
    }

    delete[] kmers;
    return false;
}


/* Returns a sorted vector by decreasing order and returns the n most 
   frequent kmers if they exist */
/*vector<int> Distance::compute_key(const string& s, int n) {
    typedef unordered_map<int,int> gmap;
    // count grams in shorter
    gmap grams;  // gram count index in lexicographical order
    int index;
    for (unsigned int i = 0; i <= s.length()-k; i++) {
        index = gram_pos(s.substr(i,k));

        if (grams.find(index) == grams.end()) {
            grams.insert(pair<int,int>(index, 1));
        }
        else {
            grams[index]++;
        }
    }

    vector<pair<int, int>> map_pairs;
    vector<int> ret;
    int i = 0;
    for (gmap::const_iterator it = grams.begin(); it != grams.end(); ++it) {
        map_pairs.push_back({it->first,it->second});
    }
    sort(map_pairs.begin(), map_pairs.end(), 
            [](const pair<int, int>& lhs, const pair<int, int>& rhs) {
                return lhs.second > rhs.second;
            });
    for (vector<pair<int, int>>::const_iterator it = map_pairs.begin(); 
            it != map_pairs.end() && i < n; ++i,++it) {
        ret.push_back(it->first);
    }
    return ret;
}*/


int Distance::levenshtein(string s, string t) {
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
    for (int i = 0; i < slen+1; i++) {
        pcol[i] = i;
    }
    // Dynamic approach to calculate the distance between two strings
    for (int i = 0; i < slen; i++) {
        col[0] = i+1;
        for (int j = 0; j < tlen; j++) {
            int cost = !(s[j] == t[i]);
            col[j+1] = min(col[j] + 1, min(pcol[j+1] + 1, pcol[j] + cost));
        }
        for (int j = 0; j < slen + 1; j++) {
            pcol[j] = col[j];
        }
    }
    return col[slen];
}

double Distance::levenshtein_window(string s, string t) {
    int slen = s.length();
    int tlen = t.length();
    
    string shorter, longer;
    int short_len, long_len, cur_dist;
    int min_dist = 9999;

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

    if (windows == 0)
        return (double)(win_size - levenshtein(s, t)) / (double)win_size;

    for (int i = 0; i < windows; i++) {
        cur_dist = levenshtein(shorter, longer.substr(i, win_size));
        min_dist = min(min_dist, cur_dist);
    }
    
    return (double)(win_size - min_dist) / (double)win_size;
}

/*void Distance::printDistMatrix(const string& filename, int k, int count, int threshold) {
    fstream fs0(filename);
    fstream fs1(filename);

    string fst, snd;
    int distances[count][count];

    for (int i = 0; i < count; i++)
        for (int j = 0; j < count; j++)
            distances[i][j] = -1; // initialize distance matrix entries to -1

    for (int i = 0; i < count; ++i) {
        IO::read_sequence(fs0, fst);

        for (int j = 0; j < count; ++j) {
            IO::read_sequence(fs1, snd);
            if (i == j) { // don't compare a sequence to itself
                distances[i][j] = 0;
                continue;
            }
            if (distances[i][j] != -1) {
                continue;
            }

            int newdist = Distance::d2window(fst, snd, k, threshold);
            distances[i][j] = newdist;
            distances[j][i] = newdist;
        }
        fs1.seekg(0, ios::beg); // rewind fs1 to start
    }

    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < count; ++j) {
            cout << setw(4) << distances[i][j]; // TODO: reset width?
        }
        cout << endl;
    }

    fs0.close();
    fs1.close();
}*/


/**
 * Calculate (0-based) index for a given k-gram in lexicographical order and
 * based on lenght k of given input string. Example:
 *   gram_pos("ag") == 2
 *   gram_pos("ca") == 4
 */
/*inline bitset<32> Distance::gram_pos(const string& s) {
    int slen = s.length();
    bitstring index;
    for (int i = 0, j = 0; i < slen; ++i, j+=2) {
        switch (s[i]) {
            case 'a':
            case 'A':
                break;
            case 'c':
            case 'C':
                index.set(j);
                break;
            case 'g':
            case 'G':
                index.set(j+1);
                break;
            case 't':
            case 'T':
            case 'u':
            case 'U':
                index.set(j);
                index.set(j+1);
                break;
            default:
                //cout << "Unknown char passed to gram_pos" << '\n';
                ;
        }
    }
    return index;
}*/

/*set<string> Distance::kmers(const Seq& s) {
    set<string> kmers;
    for (unsigned int i = 0; i <= s.data.length()-k; i++) {
        kmers.insert(s.data.substr(i,k));
    }
    return kmers;
}*/
