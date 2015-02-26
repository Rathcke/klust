#include "Distance.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
using namespace std;

int Distance::d2(string s, string t, int k) {
    int gramlen = pow(4,k);
    int slen = s.length(), tlen = t.length();
    int dist = 0;
    int v0[gramlen], v1[gramlen];
    int i;
    // Initialize two array to count how many times each substring occurs in
    // string
    for (i = 0; i < gramlen; i++) {
        v0[i] = 0;
        v1[i] = 0;
    }
    // Amount of each substring in s
    for (i = 0; i <= slen-k; i++) {
        v0[gram_pos(s.substr(i,k))]++;
    }
    // Amount of each substring in t
    for (i = 0; i <= tlen-k; i++) {
        v1[gram_pos(t.substr(i,k))]++;
    }
    // Euclidian distance between the two arrays
    for (i = 0; i < gramlen; i++) {
        dist += pow(v0[i]-v1[i],2);
    }
    return sqrt(dist);
}

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
			int cost = !(s[i] == t[j]);
			col[j+1] = min(col[j] + 1, min(pcol[j+1] + 1, pcol[j] + cost));
		}
		for (int j = 0; j < slen + 1; j++) {
			pcol[j] = col[j];
		}
	}
	return col[slen];
}
