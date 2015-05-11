#include <climits>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <list>
#include <string>
#include <vector>
#include <iostream>

#include "Cluster.h"
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

    if (count != INT_MAX)
        seqs.reserve(count); // reserve space if sequence count is specified

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
    data[seq_len] = buf[seq_len];    
    seq_len += 1;
    if (p != data) {
        seqs.emplace_back(data, seq_len, desc);
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

/**
 * Given a list of Centroid structs and two file streams, centroids and
 * clusters, write centroids in FASTA format to the former and clustering
 * results to the latter (in a format inpired by UCLUST's .uc format).
 */
void write_results(const list<Centroid>& cts, ofstream& fs_centroids,
        ofstream& fs_clusters) {

    for (auto& c : cts) {
        // write FASTA format to centroids file
        fs_centroids << '>' << c.seq.desc << '\n'
                     << c.seq.to_string() << '\n';

        // write centroid entry to clusters file
        fs_clusters << 'C' << setw(6) << c.num << ' ' << c.seq.desc << '\n';

        // write hit entries to clusters file
        for (auto& q : c.cls_seqs) {
            fs_clusters << 'H' << setw(6) << c.num
                        << ' ' << q.get().desc
                        << ' ' << c.seq.desc << '\n';
        }
    }
}

} // namespace IO
