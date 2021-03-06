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

        kmer_l = stream2int(longer + (i/4));
        kmer_l >>= (32 - 2*(i % 4) - k2);
        kmer_l &= mask;

        kmer_s = stream2int(shorter + (i/4));
        kmer_s >>= (32 - 2*(i % 4) - k2);
        kmer_s &= mask;

        kmers[kmer_l]++ < 0 ? --cur_dist : ++cur_dist;
        kmers[kmer_s]-- > 0 ? --cur_dist : ++cur_dist;

    }

    int min_dist = cur_dist;    // the least distance window so far
    int win_size = short_len;
    int windows = long_len - short_len;
    int total = 2 * (short_len - k + 1);

    double jaccard_dist = 0;

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
    }

    jaccard_dist = (double) (total - min_dist) / (double) total;

    delete[] kmers;
    return jaccard_dist;
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
