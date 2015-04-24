#ifndef IO_H
#define IO_H

#include <bitset>
#include <string>
#include <fstream>

class Seq {
    public:
        Seq();
        Seq(const char *seq_str, size_t len, const std::string& description);
        Seq(const Seq& seq);
        ~Seq();
        Seq& operator= (const Seq& seq);

        uint8_t* data() const { return seq_data; }
        size_t length() const { return seq_len; }
        //std::string& desc() { return description; }

        std::string to_string() const;

        std::string desc;

    private:
        uint8_t *seq_data;
        size_t seq_len;
        //std::string description;
        unsigned int bytes;
};

class IO
{
    public:
        IO() {}

        static void read_seqs(std::ifstream &fs, std::vector<Seq>& seqs, int count);

        static int read_seqs(std::ifstream& fs,
                std::vector<std::vector<std::bitset<2>>>& seqs, int count);

        static bool read_sequence(std::ifstream& fs, std::string& s);
        //static bool read_sequence(std::fstream& fs, struct Seq& s);

        //static void sort_incr_len(std::fstream& fs_in, std::fstream& fs_out);

    private:

};

#endif
