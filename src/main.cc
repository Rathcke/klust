#include "Distance.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

Distance dist;

int main()
{

    string a, b, c, d;
    ifstream myfile ("data.txt");
    if (myfile.is_open()) {
        getline(myfile, a);
        getline(myfile, b);
        getline(myfile, c);
        getline(myfile, d);
        myfile.close();
    } else
        return 1;

    //int d2_dist  = dist.d2(a, b, 2);            // d2 distance using 2-grams
    //int lev_dist = dist.levenshtein(a, b);      // Levenshtein distance

    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "c = " << c << endl;
    cout << "d = " << d << endl;

    cout << "d2 distance(a,b): "            << dist.d2(a,b,3)           << endl;
    cout << "Levenshtein distance(a,b): "   << dist.levenshtein(a,b)    << endl;
    cout << endl;
    cout << "d2 distance(c,d): "            << dist.d2(c,d,3)           << endl;
    cout << "Levenshtein distance(c,d): "   << dist.levenshtein(c,d)    << endl;

    return 0;
}
