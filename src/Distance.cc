#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>

#include "Distance.h"

using namespace std;

/**
 * Given two strings and an integer k, calculated the d2-distance between the
 * strings based on k-grams: count the occurences of each poosible k-mer and
 * return the Euclidean distance between the two k-mer occurence vectors.
 */
int Distance::d2(string s, string t, int k) {
    int slen = s.length(),  // length of input strings
        tlen = t.length();
    int dist = 0;           // resulting distance
    map<int, int> grams;    // gram count index in lexicographical order

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
    for (map<int,int>::iterator it = grams.begin(); it != grams.end(); ++it) {
        dist += pow(it->second, 2);
        cout << it->second << " ";
    }
    cout << endl;

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
                break;
            case 'c':
                cost += pow(4,slen-i-1);
                break;
            case 'g':
                cost += 2*pow(4,slen-i-1);
                break;
            case 't':
            case 'u':
                cost += 3*pow(4,slen-i-1);
                break;
            default:
                cout << "Unknown char passed to gram_pos" << '\n';
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
