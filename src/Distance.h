#ifndef DISTANCE_H
#define DISTANCE_H

#include <string>

class Distance
{
    public:
        Distance() {}

        int d2(std::string s, std::string t, int k);
        int levenshtein(std::string s, std::string t);
        int lev(std::string s, std::string t);

    private:
        int gram_pos(std::string s);
};

#endif
