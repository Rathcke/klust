#ifndef DISTANCE_H
#define DISTANCE_H

#include <string>

class Distance
{
    public:
        Distance() {}

        static int d2(const std::string s, const std::string t, int k);
        static bool d2window(const std::string s, const std::string t, int k, int threshold);
        static std::vector<int> computeKey(const std::string& s, int k, int n);
        static int d2window_naive(std::string s, std::string t, int k);
        static int levenshtein(std::string s, std::string t);
        static int lev(std::string s, std::string t);

        static void printDistMatrix(const std::string& filename, int k, int count, int threshold);

    private:
        static int gram_pos(std::string s);
};

#endif
