#ifndef IO_H
#define IO_H

#include <bitset>
#include <string>
#include <fstream>

struct Seq {
    std::string desc; // description following '>' in fasta format
    std::string data; // actual sequence data
};

class IO
{
    public:
        IO() {}

        static int read_seqs(std::ifstream& fs,
                std::vector<std::vector<std::bitset<2>>>& seqs, int count);

        static bool read_sequence(std::fstream& fs, std::string& s);
        static bool read_sequence(std::fstream& fs, struct Seq& s);

        static void sort_incr_len(std::fstream& fs_in, std::fstream& fs_out);

    private:

};

#endif
