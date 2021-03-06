#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Seq.h"

using namespace std;

Seq::Seq() {
    length = 0;
    bytes = 0;
    data = nullptr;
}

Seq::Seq(const char *seq_str, size_t len, const string& description) {
    length = len;
    desc = description;

    // calculate bytes necessary to store sequence as 2 bit characters in an
    // array of 8 bit unsigned ints, memset to 0 and allocate space on the heap
    bytes = len / 4 + (len % 4 != 0);
    data = new uint8_t[bytes];
    memset(data, 0, bytes);

    for (size_t i = 0; i < len; ++i) {
        int shift = 6 - 2 * (i % 4);

        switch (seq_str[i]) {
            case 'c': case 'C':
                data[i/4] |= 1 << shift; //0b01
                break;
            case 'g': case 'G':
                data[i/4] |= 2 << shift; //0b10
                break;
            case 't': case 'T': case 'u': case 'U':
                data[i/4] |= 3 << shift; //0b11
                break;
            default:
                // TODO: handle more gracefully
                break;
        }
    }
}

Seq::Seq(const Seq& seq) {
    //out << "copy cons called" << endl;
    data = new uint8_t[seq.bytes];
    memcpy(data, seq.data, seq.bytes);
    length = seq.length;
    bytes = seq.bytes;
    desc= seq.desc;
}

Seq::~Seq() {
    if (data != NULL)
        delete[] data;
}

Seq& Seq::operator= (const Seq& seq) {
    // TODO cout << "op= called" << endl;
    if (this != &seq) {
        if (data != NULL)
            delete[] data;
        data = new uint8_t[seq.bytes];
        memcpy(data, seq.data, seq.bytes);
        length = seq.length;
        bytes = seq.bytes;
        desc= seq.desc;
    }
    return *this;
}

string Seq::to_string() const {
    string s;
    s.resize(length);

    for (size_t i = 0; i < length; ++i) {
        int shift = 6 - 2 * (i % 4);
        uint8_t mask = 3 << shift;
        uint8_t bit = (data[i/4] & mask) >> shift;

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

/**
 * TODO: this might be pretty cool and useful!
 *
 * Extract substring, i.e k-mer of k = len, from the given position in the
 * current sequence object and of the given length.
 */
uint32_t Seq::substr(size_t pos, size_t len) const {
    static const int kmer_count = pow(4, len);
    static const uint32_t k2 = 2 * len;
    static const uint32_t mask = kmer_count - 1; // 0b001111 (2*k 1's)

    uint32_t kmer = 0;  // binary repr. of substring in sequence
    size_t index = pos / 4; // position in uint8_t array

    kmer = (uint32_t) data[index]   << 24 |
           (uint32_t) data[index+1] << 16 |
           (uint32_t) data[index+2] <<  8 |
           (uint32_t) data[index+3];

    // right shift k-mer to least significant byte (i.e. lowest possible value)
    kmer >>= (32 - 2*(pos % 4) - k2);
    kmer &= mask;

    return kmer;
}
