#include <climits>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <list>
#include <string>
#include <vector>

#include "Cluster.h"
#include "IO.h"
#include "Seq.h"

namespace IO {

using namespace std;

int read_seqs2(istream &fs, vector<Seq>& seqs, int count) {
    if (count != INT_MAX)
        seqs.reserve(count); // reserve space if sequence count is specified

    string seq, desc;

    int i = 0;
    while (i < count) {
        seq.clear();
        desc.clear();

        if (fs)
            getline(fs, desc);
        else
            break;
        getline(fs, seq, '>');

        seqs.emplace_back(seq.c_str(), seq.length(), desc);
        ++i;
    }

    return seqs.size();
}

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
    if (p != data) {
        seqs.emplace_back(data, seq_len, desc);
        p = data;
        seq_len = 0;
        i++;
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

void springy(const list<Centroid>& cts, ofstream& fs) {

    fs << "<html>\n"
          "<body>\n"
          "<script src=\"http://ajax.googleapis.com/ajax/libs/jquery/"
          "1.3.2/jquery.min.js\"></script>\n"
          "<script src=\"springy.js\"></script>\n"
          "<script src=\"springyui.js\"></script>\n"
          "<script>\n"
          "var graph = new Springy.Graph();\n\n";

    int cent_num = 0;
    int seq_num = 0;

    for (auto it = cts.cbegin(); it != cts.cend(); ++it) {
        fs << "var c" << cent_num
           << " = graph.newNode({label: 'c" << cent_num << "'});\n";

        vector<reference_wrapper<const Seq>> seqs = it->cls_seqs;
        for (auto s_it = seqs.cbegin(); s_it != seqs.cend(); ++s_it) {
            fs << "var s" << seq_num
               << " = graph.newNode({label: 's" << seq_num << "'});\n"
               << "graph.newEdge(c" << cent_num << ", s" << seq_num << ", "
                                    << "{color: '#0000ff'});\n";
            ++seq_num;
        }
        ++cent_num;
    }

    fs << "\njQuery(function(){\n"
       << "  var springy = window.springy = jQuery('#springydemo').springy({\n"
       << "    graph: graph,\n"
       << "    nodeSelected: function(node){\n"
       << "      console.log('Node selected: ' + JSON.stringify(node.data));\n"
       << "    }\n"
       << "  });\n"
       << "});\n"
       << "</script>\n\n"
       << "<canvas id=\"springydemo\" width=\"1250\" height=\"650\" />\n"
       << "</body>\n"
       << "</html>\n";
}

} // namespace IO
