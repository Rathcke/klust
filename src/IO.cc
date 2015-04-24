#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "IO.h"

using namespace std;

Seq::Seq() {
    seq_len = 0;
    bytes = 0;
    seq_data = NULL;
}

Seq::Seq(const char *seq_str, size_t len, const string& description) {
    seq_len = len;
    desc = description;

    // calculate bytes necessary to store sequence as 2 bit characters in an
    // array of 8 bit unsigned ints, memset to 0 and allocate space on the heap
    bytes = len / 4 + (len % 4 != 0);
    seq_data = new uint8_t[bytes];
    memset(seq_data, 0, bytes);

    for (size_t i = 0; i < len; ++i) {
        int shift = 6 - 2 * (i % 4);

        switch (seq_str[i]) {
            case 'c': case 'C':
                seq_data[i/4] |= 1 << shift; //0b01
                break;
            case 'g': case 'G':
                seq_data[i/4] |= 2 << shift; //0b10
                break;
            case 't': case 'T': case 'u': case 'U':
                seq_data[i/4] |= 3 << shift; //0b11
                break;
            default:
                // TODO: handle more gracefully
                break;
        }
    }
    //cout << bitset<8>(seq_data[0]) << endl;;
}

Seq::Seq(const Seq& seq) {
    //out << "copy cons called" << endl;
    seq_data = new uint8_t[seq.bytes];
    memcpy(seq_data, seq.seq_data, seq.bytes);
    seq_len = seq.seq_len;
    bytes = seq.bytes;
    desc = seq.desc;
}

Seq::~Seq() {
    if (seq_data != NULL)
        delete[] seq_data;
}

Seq& Seq::operator= (const Seq& seq){
    // TODO cout << "op= called" << endl;
    if (this != &seq) {
        if (seq_data != NULL)
            delete[] seq_data;
        seq_data = new uint8_t[seq.bytes];
        memcpy(seq_data, seq.seq_data, seq.bytes);
        seq_len = seq.seq_len;
        bytes = seq.bytes;
        desc= seq.desc;
    }
    return *this;
}

string Seq::to_string() const {
    string s;
    s.resize(seq_len);

    for (size_t i = 0; i < seq_len; ++i) {
        int shift = 6 - 2 * (i % 4);
        uint8_t mask = 3 << shift;
        uint8_t bit = (seq_data[i/4] & mask) >> shift;

        switch (bit) {
            case 0:
                s[i] = 'A';
                break;
            case 1:
                s[i] = 'C';
                break;
            case 2:
                s[i] = 'G';
                break;
            case 3:
                s[i] = 'T';
                break;
            default:
                throw invalid_argument("invalid bit representation "
                                       "in compressed sequence");
                break;
        }
    }
    return s;
}


void IO::read_seqs(ifstream &fs, vector<Seq>& seqs, int count) {
    // replacing internal buffer with a larger one, might speed things up
    //const unsigned int in_buf_size = 1024 * 1024;
    //char *in_buf = new char[in_buf_size];
    //fs_in.rdbuf()->pubsetbuf(in_buf, in_buf_size);
    // (...)
    //delete[] in_buf;

    seqs.reserve(count);

    string desc;

    // sequence data buffer
    const size_t data_size = 16 * 1024;
    char *data = new char[data_size];
    char *p = data;
    size_t seq_len = 0;

    // reading buffer
    const size_t buf_size = 16 * 1024;
    char *buf = new char[buf_size];

    streamsize bytes_read;

    int i = 0;
    while (fs.getline(buf, buf_size) && i < count) {
        if (buf[0] == '>') {
            if (p != data) {
                seqs.emplace_back(data, seq_len, desc);
                p = data;
                seq_len = 0;
                i++;
            }
            buf[fs.gcount() - 1] = '\0';
            desc.assign(buf + 1);
            continue;
        }

        // TODO: fix possibility for buffer overflow
        bytes_read = fs.gcount() - 1; // gcount() counts discarded \n as well
        memcpy(p, buf, bytes_read);
        p += bytes_read;
        seq_len += bytes_read;
    }
}

/**
 * Read [DR]NA sequence from given filestream in given string.
 * Return true on success and false if there's no more to read.
 */
bool IO::read_sequence(ifstream& fs, string& s) {
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
/*void sort_incr_len(fstream& fs_in, fstream& fs_out) {
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
}*/
