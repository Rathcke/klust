#ifndef DISTANCE_H
#define DISTANCE_H

#include <string>
#include <set>
#include <bitset>
#include "IO.h"

class Distance
{
    public:
        Distance(int kmer, double threshold, int step_size);

        bool compare(const Seq& s, const Seq& t);

        bool compare(const std::vector<std::bitset<2>>& s,
                const std::vector<std::bitset<2>>& t);

        bool compare(const std::string& s, const std::string& t);

        std::vector<int> compute_key(const std::string& s, int n);

        double levenshtein(const std::string& s, const std::string& t);

        double levenshtein_window(const std::string& s, const std::string& t);

        std::set<std::string> kmers(const Seq& s);

        void printDistMatrix(const char* filename, int count);

    private:
        int k;          // k in k-mer (word length)
        double thrs;    // threshold
        int step;
};

#endif
