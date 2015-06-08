#ifndef DISTANCE_H
#define DISTANCE_H

#include <string>
#include <set>

#include "Seq.h"

class Distance
{
    public:
        Distance(int kmer, double threshold)
            : k {kmer}, thrs {threshold} {}

        /**
         * Compare two Seqs and return true if their similarity is greater than
         * or equal to the threshold similarity.
         */
        bool compare(const Seq& s, const Seq& t);

        /**
         * Calculate the k-mer based similarity between the two given Seqs.
         * Return a value in the interval [0,1].
         */
        double distance(const Seq& s, const Seq& t);

        /**
         * Returns the Levenshtein distance between the two given strings.
         * Uses a bottom-up dynamic programming algorithm.
         */
        double levenshtein(const std::string& s, const std::string& t);

        /**
         * Return the threshold.
         */
        inline double threshold() { return thrs; }

        /**
         * Return the value for the distance parameter k.
         */
        inline int kmer() { return k; }

        /**
         * Given a pointer to an array of uint8_t, extract a uint32_t
         * containing the bytes from four uint8_t in a row.
         */
        static inline uint32_t stream2int(const uint8_t *stream) {
            return (((uint32_t) stream[0]) << 24 |
                    ((uint32_t) stream[1]) << 16 |
                    ((uint32_t) stream[2]) <<  8 |
                    ((uint32_t) stream[3]) <<  0);
        }

    private:
        int k;          // k in k-mer (word length)
        double thrs;    // threshold
};

#endif
