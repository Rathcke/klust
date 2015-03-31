#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "IO.h"
#include "Distance.h"
#include "Cluster.h"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 9) {
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
    }

    fstream fs_in(argv[1]);
    /*fstream fs_out(argv[2], fstream::out | fstream::trunc);
    fstream fs_cts(argv[2], fstream::out | fstream::trunc);
    fstream fs_cls(argv[3], fstream::out | fstream::trunc);
    const int k = std::atoi(argv[4]);         // k in k-mers
    const double threshold = stod(argv[5]); // simlilarity threshold
    const int count = std::atoi(argv[6]);     // # of sequences to measure
    const int max_rejects = std::atoi(argv[7]);
    const int step_size = std::atoi(argv[8]);

    Distance d2(k, threshold, step_size);*/


    int count = 1;
    vector<vector<bitset<2>>> seqs;
    clock_t start = clock();
    IO::read_seqs(fs_in, seqs, count);
    double a = (clock()-start) / (double)(CLOCKS_PER_SEC);

    for (unsigned int i = 0; i < seqs.size(); ++i) {
        for (unsigned int j = 0; j < seqs[i].size(); ++j) {
            cout << seqs[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Time: " << a << endl;
    cout << "Seqs/sec: " << count / a << endl;

    /*Seq s;
    vector<Seq> seqs;
    for (int i = 0; i < count; ++i) {
        IO::read_sequence(fs_in, s);
        seqs.push_back(s);
    }

    clock_t start = clock();
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < count; ++j) {
            d2.compare(seqs[i].data, seqs[j].data);
        }
    }
    double a = (clock()-start) / (double)(CLOCKS_PER_SEC);
    cout << "Time: " << a << endl;
    cout << "Comparisons/sec: " << pow(count, 2) / a << endl;
    return 0;

    cout << "# of clusters: " <<
        //Cluster::clust(fs_in, fs_cts, threshold, k, count) << endl;
        Cluster::intersect_clust(fs_in, fs_cts, fs_cls, d2, count, max_rejects) << endl;

    //Distance::printDistMatrix(argv[1], k, count);

    fs_in.close();
    fs_out.close();
    fs_cts.close();
    fs_cls.close();*/

    return 0;
}
