#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <fstream>

#include "IO.h"

using namespace std;

/**
 * Read [DR]NA sequence from given filestream in given string.
 * Return 0 on success and -1 if there's no more to read.
 */
bool IO::readSequence(fstream& fs, string& s) {
    string tmp;
    s = ""; // clear s

    getline(fs, tmp);
    if (tmp.empty())
        return false;
    while (tmp[0] == '>')
        getline(fs, tmp);                  // ignore '>' line
    while (!tmp.empty() && tmp[0] != '>') {     // read until next '>' line
        s += tmp;
        getline(fs, tmp);
    }
    return true;
}