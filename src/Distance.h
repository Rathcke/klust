#ifndef DISTANCE_H
#define DISTANCE_H

#include <string>

using namespace std;

class Distance
{
    public:
        Distance() {}

        int d2(string s, string t, int k);
        int levenshtein(string s, string t);

    private:
        int gram_pos(string s);
};

#endif
