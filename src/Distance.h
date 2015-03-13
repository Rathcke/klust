#ifndef DISTANCE_H
#define DISTANCE_H

#include <string>

class Distance
{
    public:
        Distance(int kmer, double threshold);

        bool compare(const std::string& s, const std::string& t);

        static std::vector<int> computeKey(const std::string& s, int k, int n);

        static int levenshtein(std::string s, std::string t);

        /*static void printDistMatrix(const std::string& filename,
                int k, int count, int threshold);*/

    private:
        int k;      // k in k-mer (word length)
        double thrs;   // threshold

        static int gram_pos(const std::string& s);


};

#endif
