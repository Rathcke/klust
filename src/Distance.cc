#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "Distance.h"
#include "IO.h"

using namespace std;

Distance::Distance(int kmer, double threshold) {
    this->k = kmer;
    this->thrs = threshold;
}

/**
 * Given two strings, the int k in k-mer (word length), calculated the
 * d2-distance between the strings based on k-grams: count the occurences of
 * each poosible k-mer, calculate the Euclidean distance between the two k-mer
 * occurence vectors. Return true if within threshold in some window (of size
 * shortest string), otherwise return false.
 */
bool Distance::compare(const string& s, const string& t) {
    int slen = s.length(),
        tlen = t.length();
    string shorter, longer;
    int short_len, long_len;
    int total = 0;

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

    typedef unordered_map<int,int> gmap;
    // count grams in shorter
    gmap grams; // gram count index in lexicographical order

    int index;
    for (int i = 0; i <= short_len-k; i++) {
        index = gram_pos(shorter.substr(i,k));

        if (grams.find(index) == grams.end()) {
            grams.insert(pair<int,int>(index, 1));
        }
        else {
            grams[index]++;
        }
        total += 1;
    }


    // Amount of each substring in t
    for (int i = 0; i <= short_len-k; i++) {
        index = gram_pos(longer.substr(i,k));

        if (grams.find(index) == grams.end()) {
            grams.insert(pair<int,int>(index, -1));
        }
        else {
            grams[index]--;
        }
        total += 1;
    }

    int init = 0;
    // euclidian distance between the two arrays
    for (gmap::iterator it = grams.begin(); it != grams.end(); ++it)
        init += abs(it->second);

    int min_dist = init; // variable containing the least distance window so far
    int cur_dist = init; // distance in current window

    int win_size = short_len;
    int windows = long_len - short_len;

    if (windows == 0)
        return cur_dist >= thrs;

    string pre_gram, post_gram;
    for (int i = 0; i < windows; i++) {
        pre_gram  = longer.substr(i, k);
        post_gram = longer.substr(i + win_size - k + 1, k);

        if (pre_gram == post_gram)
            continue;

        // If changed for the better, decrement cur_dist, otherwise increment it.
        grams[gram_pos(pre_gram)] < 0 ? --cur_dist : ++cur_dist;
        grams[gram_pos(post_gram)] > 0 ? --cur_dist : ++cur_dist;

        grams[gram_pos(pre_gram)] += 1;

        int post_gram_pos = gram_pos(post_gram);
        if (grams.find(post_gram_pos) == grams.end())
            grams.insert(pair<int,int>(post_gram_pos, -1));
        else
            grams[gram_pos(post_gram)] -= 1;

        min_dist = min(min_dist, cur_dist);

        if ((double)(total-min_dist)/(double)total >= thrs) {
            // cout << (double)(total-min_dist)/(double)total << endl;
            return true;
        }
    }

    return false;
}

/* Returns a sorted vector by decreasing order and returns the n most 
   frequent kmers if they exist */
vector<pair<int,int>> Distance::compute_key(const string& s, unsigned int n) {
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
    //vector<int> ret;
    for (gmap::const_iterator it = grams.begin(); it != grams.end(); ++it) {
        map_pairs.push_back({it->first,it->second});
    }
    sort(map_pairs.begin(), map_pairs.end(), 
            [](const pair<int, int>& lhs, const pair<int, int>& rhs) {
                return lhs.second > rhs.second;
            });
    if (map_pairs.size() >= n)
        map_pairs.resize(n);
    return map_pairs;

/*    for (vector<pair<int, int>>::const_iterator it = map_pairs.begin(); 
            it != map_pairs.end() && i < n; ++i,++it) {
        ret.push_back(it->first, it->second);
    }
    return ret;*/
}


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
    for (int i = 0; i < tlen; i++) {
        col[0] = i+1;
        for (int j = 0; j < slen; j++) {
            int cost = !(s[j] == t[i]);
            col[j+1] = min(col[j] + 1, min(pcol[j+1] + 1, pcol[j] + cost));
        }
        for (int j = 0; j < slen + 1; j++) {
            pcol[j] = col[j];
        }
    }
    return col[slen];
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
int Distance::gram_pos(const string& s) {
    int slen = s.length();
    int cost = 0;
    // Loop that calculates the index for a substring
    for (int i = slen - 1; i >= 0; i--) {
        switch (s[i]) {
            case 'a':
            case 'A':
                break;
            case 'c':
            case 'C':
                cost += pow(4,slen-i-1);
                break;
            case 'g':
            case 'G':
                cost += 2*pow(4,slen-i-1);
                break;
            case 't':
            case 'T':
            case 'u':
            case 'U':
                cost += 3*pow(4,slen-i-1);
                break;
            default:
                //cout << "Unknown char passed to gram_pos" << '\n';
                ;
        }
    }
    return cost;
}
