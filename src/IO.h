#ifndef IO_H
#define IO_H

#include <bitset>
#include <string>
#include <fstream>

/*struct Seq {
    std::string desc; // description following '>' in fasta format
    std::string data; // actual sequence data
};*/

class Seq {
    public:
        Seq(const char *seq_str, size_t len);
        Seq(const Seq& seq);
        ~Seq();
        Seq& operator= (const Seq& seq);

        uint8_t* data() const { return seq_data; }
        size_t length() const { return seq_len; }

        std::string to_string() const;

    private:
        uint8_t *seq_data;
        size_t seq_len;
        unsigned int bytes;
};

class IO
{
    public:
        IO() {}

        void read_seqs(std::ifstream &fs, std::vector<Seq>& seqs, int count);

        static int read_seqs(std::ifstream& fs,
                std::vector<std::vector<std::bitset<2>>>& seqs, int count);

        static bool read_sequence(std::fstream& fs, std::string& s);
        //static bool read_sequence(std::fstream& fs, struct Seq& s);

        //static void sort_incr_len(std::fstream& fs_in, std::fstream& fs_out);

    private:

};

#endif
