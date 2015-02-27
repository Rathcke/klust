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
    cout << endl;

    //string a = "aaaaagg";
    //string b = "cccccgg";
    string e = "abcbadfe";
    string f = "bacbdfe";

    cout << "d2 distance(a,b): "            << dist.d2(a,b,2)           << endl;
    cout << "Levenshtein distance(a,b): "   << dist.levenshtein(a,b)    << endl;

    //string a = "acgtgtgacgtgcatagcgtacgtgac";
    //string b = "gtgtacgatagtcgtactgatcatg";

    return 0;
}
