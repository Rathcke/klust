#ifndef SEQ_H
#define SEQ_H

#include <bitset>

#define KMER_BITSET 4096
#define KMER_LEN 6

class Seq {
    public:
        Seq();
        Seq(const char *seq_str, size_t len, const std::string& description);
        Seq(const Seq& seq);
        ~Seq();
        Seq& operator= (const Seq& seq);

        inline uint8_t* data() const { return seq_data; }
        inline size_t length() const { return seq_len; }
        const std::string& desc() const { return description; }

        std::string to_string() const;

        std::bitset<KMER_BITSET> bits;  // bitset of the k-mers occurring in seq
        size_t count;                   // # of distinct k-mers in seq; 

        void set_bits(std::bitset<KMER_BITSET>& b) {
            bits = b;
            count = b.count();
        }

        //std::string desc;

        uint32_t substr(size_t pos, size_t len) const;

    private:
        uint8_t *seq_data;
        size_t seq_len;
        std::string description;
        unsigned int bytes;

};

#endif
