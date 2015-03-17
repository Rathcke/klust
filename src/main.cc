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
    if (argc < 7) {
        std::cout << "Usage: " << argv[0] << " <.fasta input file> "
                                             " <.fasta output file for centroids> "
                                             " <output file for clusters> "
                                             " <k in k-mers> "
                                             " <similarity threshold>"
                                             " <# of sequences to compare>"
                                          << endl << endl;
        return 1;
    }

    fstream fs_in(argv[1]);
    fstream fs_out(argv[2], fstream::out | fstream::trunc);
    fstream fs_cts(argv[2], fstream::out | fstream::trunc);
    fstream fs_cls(argv[3], fstream::out | fstream::trunc);
    const int k = std::atoi(argv[4]);         // k in k-mers
    const double threshold = stod(argv[5]); // simlilarity threshold
    const int count = std::atoi(argv[6]);     // # of sequences to measure

    Distance d2(k, threshold);

    cout << "# of clusters: " <<
        //Cluster::clust(fs_in, fs_cts, threshold, k, count) << endl;
        Cluster::clust(fs_in, fs_cts, fs_cls, d2, count) << endl;

    //Distance::printDistMatrix(argv[1], k, count);

    fs_in.close();
    fs_out.close();
    fs_cts.close();
    fs_cls.close();

    return 0;
}
