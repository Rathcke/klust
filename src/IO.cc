#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <fstream>

#include "IO.h"

using namespace std;

/**
 * Given a char, return a bitset representation of the char.
 */
inline bitset<2> gram_to_bitset(const char c) {
    switch(c) {
        case 'a': case 'A':
            return bitset<2>();  // 0b00
            break;
        case 'c': case 'C':
            return bitset<2>(1); // 0b01
            break;
        case 'g': case 'G':
            return bitset<2>(2); // 0b10
            break;
        case 't': case 'T': case 'u': case 'U':
            return bitset<2>(3); // 0b11
            break;
        default:
            //cout << "unknown char passed to gram_to_bitset: " << c << endl;
            return bitset<2>();  // 0b00
            break;
    }
}

/**
 * Read specified number of sequences from filestream and load into given
 * vector. Discard the '>' lines. Return the number of sequences read.
 */
int IO::read_seqs(ifstream& fs, vector<vector<bitset<2>>>& seqs, int count) {
    const unsigned int buf_size = 16 * 1024;
    char* buf = new char[buf_size];

    vector<bitset<2>> v;
    seqs.resize(count); // segfault without for some reason

    int i = 0;
    while(fs.good() && i++ < count) {
        fs.ignore(256, '\n');   // discard first 256 chars or until newline
        // read until '>', i.e. the sequence data; note this also reads \n chars
        if (!fs.getline(buf, buf_size, '>'))
            break;

        for (char *p = buf; *p != '\0'; ++p) {  // use transform?
            if (*p == '\n')
                continue;
            v.push_back(gram_to_bitset(*p));
        }

        seqs.push_back(v);
        v.clear();
    }

    delete[] buf;

    return seqs.size();
}


/**
 * Read one entry of description and sequence data from FASTA format stream,
 * place in given seq struct; return true on success and false on failure.
 */
bool IO::read_sequence(fstream& fs, struct Seq& s) {
    string tmp;

    if (fs) { // true if stream is ready; false on EOF or error
        s.desc.clear(); // clear contents of given struct argument
        s.data.clear();

        getline(fs, tmp);
        if (!tmp.empty() && tmp[0] == '>') {
            int i = tmp.find_first_not_of("> "); // remove '>' and spaces
            s.desc = tmp.substr(i, s.desc.size() - i);
        }

        while (fs && fs.peek() != '>' && fs.peek() != EOF) {
            getline(fs, tmp);
            s.data += tmp;
        }
        return true;
    }
    return false;
}


/**
 * Read [DR]NA sequence from given filestream in given string.
 * Return true on success and false if there's no more to read.
 */
bool IO::read_sequence(fstream& fs, string& s) {
    string tmp;
    s.clear();

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

/**
 * Read seqeunces from input stream, sort in increasing (non-decreasing) order
 * and output to output stream. Sorts everything in-place in memory.
 */
void sort_incr_len(fstream& fs_in, fstream& fs_out) {
    Seq s;
    vector<Seq> seqs;
    while (IO::read_sequence(fs_in, s)) {
        seqs.push_back(s);
    }

    sort(seqs.begin(), seqs.end(),
            [](Seq s1, Seq s2) {
                return s1.data.length() < s2.data.length();
            });

    for(vector<Seq>::const_iterator it = seqs.begin(); it != seqs.end(); ++it) {
        fs_out << '>' << (*it).desc << endl
               << (*it).data << endl;
    }
}
