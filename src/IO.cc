#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "IO.h"
#include "Seq.h"

namespace IO {

using namespace std;

/**
 * Read the specified number of sequences from the given filestream (in FASTA
 * format), construct a Seq object for each sequence and store in the given
 * vector. Return the number of read sequences.
 */
int read_seqs(ifstream &fs, vector<Seq>& seqs, int count) {
    // replacing internal buffer with a larger one, might speed things up
    //const unsigned int in_buf_size = 1024 * 1024;
    //char *in_buf = new char[in_buf_size];
    //fs_in.rdbuf()->pubsetbuf(in_buf, in_buf_size);
    // (...)
    //delete[] in_buf;

    //seqs.reserve(count); TODO: maybe reserve if we know count

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

    delete[] data;
    delete[] buf;

    return seqs.size();
}

/**
 * Read [DR]NA sequence from given filestream in given string.
 * Return true on success and false if there's no more to read.
 */
bool read_sequence(ifstream& fs, string& s) {
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

} // namespace IO
