#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cstring>

#include "Cluster.h"
#include "Distance.h"
#include "IO.h"

using namespace std;

int main(int argc, char *argv[])
{
    /*if (argc < 9) {
        std::cout << "Usage: " << argv[0] << " <.fasta input file> "
                                             " <.fasta output file for centroids> "
                                             " <output file for clusters> "
                                             " <k in k-mers> "
                                             " <similarity threshold>"
                                             " <# of sequences to compare>"
                                             " <max_rejects>"
                                             " <step_size>"
                                          << endl << endl;
        return 1;
    }*/

    /*cout << sizeof(bitset<65>) << endl;
    return 0;*/

    ios_base::sync_with_stdio(false); // don't share buffers with C style IO

    ifstream fs_in(argv[1]);
    if (!fs_in) {
        cerr << "error opening file: " << argv[1] << endl;
        return 1;
    }

    // replacing internal buffer with a larger one, might speed things up
    //const unsigned int in_buf_size = 1024 * 1024;
    //char *in_buf = new char[in_buf_size];
    //fs_in.rdbuf()->pubsetbuf(in_buf, in_buf_size);
    // (...)
    //delete[] in_buf;

    /*fstream fs_out(argv[2], fstream::out | fstream::trunc);
    fstream fs_cts(argv[2], fstream::out | fstream::trunc);
    fstream fs_cls(argv[3], fstream::out | fstream::trunc);
    const int k = std::atoi(argv[4]);         // k in k-mers
    const double threshold = stod(argv[5]); // simlilarity threshold
    const int count = std::atoi(argv[6]);     // # of sequences to measure
    const int max_rejects = std::atoi(argv[7]);
    const int step_size = std::atoi(argv[8]);

    Distance d2(k, threshold, step_size);*/
    Distance d2(8, 0.8, 0);

    /*char a[] = "acatgatgcagt";
    Seq s(a, strlen(a));
    Seq t(a, strlen(a));
    d2.compare(s, t);
    return 0;*/


    /*
     * Reading sequences
     */
    int count = 500;
    vector<vector<bitset<2>>> seqs;
    //vector<Seq> seqs;

    cout << "Reading " << count << " sequences...\n" << endl;
    clock_t read_clock = clock();
    IO::read_seqs(fs_in, seqs, count);
    double read_secs = (clock() - read_clock) / (double) CLOCKS_PER_SEC;

    cout << "Finished reading:\n"
         << "Time: "     << read_secs << " sec.\n"
         << "Seqs/sec: " << count / read_secs << "\n" << endl;


    /*
     * Comparing sequences
     */
    int tot = 0;
    cout << "Comparing all read sequences...\n" << endl;
    clock_t comp_clock = clock();
    for (int i = 0; i < count; ++i)
        for (int j = 0; j < count; ++j) {
            if (d2.compare(seqs[i], seqs[j]))
            	++tot;
        }    	
    double comp_secs = (clock() - comp_clock) / (double) CLOCKS_PER_SEC;

    cout << "Finished comparing:\n"
         << "Time: "            << comp_secs << " sec.\n"
         << "Comparisons/sec: " << pow(count, 2) / comp_secs << "\n"
         << "# of compares: "   << count * count << endl;
    cout << tot << endl;

    /*cout << "# of clusters: " <<
        //Cluster::clust(fs_in, fs_cts, threshold, k, count) << endl;
        Cluster::intersect_clust(fs_in, fs_cts, fs_cls, d2, count, max_rejects) << endl;

    //Distance::printDistMatrix(argv[1], k, count);

    fs_out.close();
    fs_cts.close();
    fs_cls.close();*/

    fs_in.close();

    return 0;
}
