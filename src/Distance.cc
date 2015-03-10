#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>

#include "Distance.h"
#include "IO.h"

using namespace std;

int Distance::d2window(const string s, const string t, int k) {
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
    gmap grams;  // gram count index in lexicographical order
    int index;
    for (int i = 0; i <= short_len-k; i++) {
        index = gram_pos(shorter.substr(i,k));

        if (grams.find(index) == grams.end()) {
            grams.insert(pair<int,int>(index, 1));
        }
        else {
            grams[index]++;
        }
        total += 2*grams[index]-1;  // sum of squared gram counts in shorter
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
    }

    int init = 0;
    // euclidian distance between the two arrays
    for (gmap::iterator it = grams.begin(); it != grams.end(); ++it)
        init += pow(it->second, 2);

    int min_dist = init; // variable containing the least distance window so far
    int cur_dist = init; // distance in current window

    int win_size = short_len;
    int windows = long_len - short_len;
    string pre_gram, post_gram;
    for (int i = 0; i < windows; i++) {
        pre_gram  = longer.substr(i, k);
        post_gram = longer.substr(i + win_size - k + 1, k);

        if (pre_gram == post_gram)
            continue;

        grams[gram_pos(pre_gram)] += 1;

        int post_gram_pos = gram_pos(post_gram);
        if (grams.find(post_gram_pos) == grams.end())
            grams.insert(pair<int,int>(post_gram_pos, -1));
        else
            grams[gram_pos(post_gram)] -= 1;

        cur_dist = cur_dist + 2*grams[gram_pos(pre_gram)] -
                    2*grams[gram_pos(post_gram)] - 2;

        min_dist = min(min_dist, cur_dist);
    }

    return sqrt(min_dist);
}

int Distance::d2window_naive(string s, string t, int k) {
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

    int cur_dist, min_dist = 999999; // ugly I know
    int win_size = short_len;
    int windows = long_len - short_len;

    for (int i = 0; i < windows; i++) {
        cur_dist = d2(shorter, longer.substr(i, win_size), k);
        min_dist = min(min_dist, cur_dist);
    }
    return min_dist;
}


/**
 * Given two strings and an integer k, calculated the d2-distance between the
 * strings based on k-grams: count the occurences of each poosible k-mer and
 * return the Euclidean distance between the two k-mer occurence vectors.
 */
int Distance::d2(const string s, const string t, int k) {
    int slen = s.length(),  // length of input strings
        tlen = t.length();
    int dist = 0;           // resulting distance

    typedef unordered_map<int,int> gmap;
    gmap grams;    // gram count index in lexicographical order

    int i, index;
    // Amount of each substring in s
    for (i = 0; i <= slen-k; i++) {
        index = gram_pos(s.substr(i,k));

        if (grams.find(index) == grams.end())
            grams.insert(pair<int,int>(index, 1));
        else
            grams[index]++;
    }

    // Amount of each substring in t
    for (i = 0; i <= tlen-k; i++) {
        index = gram_pos(t.substr(i,k));

        if (grams.find(index) == grams.end())
            grams.insert(pair<int,int>(index, -1));
        else
            grams[index]--;
    }

    // Euclidian distance between the two arrays
    for (gmap::iterator it = grams.begin(); it != grams.end(); ++it) {
        dist += pow(it->second, 2);
    }

    return sqrt(dist);
}

/**
 * Calculate (0-based) index for a given k-gram in lexicographical order and
 * based on lenght k of given input string. Example:
 *   gram_pos("ag") == 2
 *   gram_pos("ca") == 4
 */
int Distance::gram_pos(string s) {
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

void Distance::printDistMatrix(const string& filename, int k, int count) {
    fstream fs0(filename);
    fstream fs1(filename);

    string fst, snd;
    int distances[count][count];

    for (int i = 0; i < count; i++)
        for (int j = 0; j < count; j++)
            distances[i][j] = -1; // initialize distance matrix entries to -1

    for (int i = 0; i < count; ++i) {
        IO::readSequence(fs0, fst);

        for (int j = 0; j < count; ++j) {
            IO::readSequence(fs1, snd);
            if (i == j) { // don't compare a sequence to itself
                distances[i][j] = 0;
                continue;
            }
            if (distances[i][j] != -1) {
                continue;
            }

            int newdist = Distance::d2window(fst, snd, k);
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
}
