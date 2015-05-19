#ifndef SEQ_H
#define SEQ_H

#include <fstream>
#include <vector>

extern unsigned long long dirty;

struct Seq {
    Seq();
    Seq(const char *seq_str, size_t len, const std::string& description);
    Seq(const Seq& seq);
    ~Seq();
    Seq& operator= (const Seq& seq);

    uint8_t *data;
    size_t length;
    size_t bytes;   // # of uint8_t required for storing the data
    std::string desc;

    std::string to_string() const;

    uint32_t substr(size_t pos, size_t len) const;

    static void permute(std::vector<Seq>& seqs, int count, double ratio,
        std::ofstream& fs_cts);
};

#endif
