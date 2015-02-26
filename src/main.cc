#include "Distance.h"
#include <string>
#include <iostream>

Distance dist;

int main()
{
    std::string a = "actacgatgctagctggtcgatcgatgctag";
    std::string b = "acgtcgtccctcgatcgctacgtggctaggagatc";


    int d2_dist  = dist.d2(a, b, 2);            // d2 distance using 2-grams
    int lev_dist = dist.levenshtein(a, b);      // Levenshtein distance

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "d2 distance: " << d2_dist << std::endl;
    std::cout << "Levenshtein distance: " << lev_dist << std::endl;

    return 0;
}
