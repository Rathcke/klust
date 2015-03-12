#ifndef DISTANCE_H
#define DISTANCE_H

#include <string>

class Distance
{
    public:
        Distance(int kmer, int threshold);

        bool compare(const std::string& s, const std::string& t);

        static int d2window(const std::string& s, const std::string& t, int k);

    private:
        int k,      // k in k-mer (word length)
            thrs;   // threshold

        static int gram_pos(const std::string& s);
};

#endif
