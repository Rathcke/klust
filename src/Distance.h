#ifndef DISTANCE_H
#define DISTANCE_H

#include <string>
#include <set>

#include "Seq.h"

class Distance
{
    public:
        Distance(int kmer, double threshold, int step_size);

        bool compare(const Seq& s, const Seq& t);
        double distance(const Seq& s, const Seq& t);

        /**
         * Given a sequence, return the int representations of up to the n most
         * frequent k-mers, if they exist.
         */
        std::vector<int> compute_key(const Seq& s, int n);

        double levenshtein(const std::string& s, const std::string& t);

        double levenshtein_window(const std::string& s, const std::string& t);

        std::set<std::string> kmers(const Seq& s);

        void printDistMatrix(const char* filename, int count);

        void jac_printDistMatrix(const char* filename, int count);

        inline double threshold() { return thrs; }

        static inline uint32_t stream2int(const uint8_t *stream) {
            return (((uint32_t) stream[0]) << 24 |
                    ((uint32_t) stream[1]) << 16 |
                    ((uint32_t) stream[2]) <<  8 |
                    ((uint32_t) stream[3]) <<  0);
        }

    private:
        int k;          // k in k-mer (word length)
        double thrs;    // threshold
        int step;
};

#endif
