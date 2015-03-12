#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <assert.h>

#include "IO.h"
#include "Distance.h"
#include "Cluster.h"

using namespace std;

bool cmpSeq(Seq s1, Seq s2) {
    return s1.data.length() < s2.data.length();
}

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
    const int threshold = std::atoi(argv[5]); // simlilarity threshold
    const int count = std::atoi(argv[6]);     // # of sequences to measure

    cout << "# of clusters: " <<
        //Cluster::clust(fs_in, fs_cts, threshold, k, count) << endl;
        Cluster::clust(fs_in, fs_cts, fs_cls, threshold, k, count) << endl;

    //Distance::printDistMatrix(argv[1], k, count);

    /* SORTING BY INCREASING LENGTH:
    fstream fs_in(argv[1]);
    fstream fs_out(argv[2], fstream::out | fstream::trunc);

    Seq s;
    vector<Seq> seqs;
    while (IO::readSequence(fs_in, s)) {
        seqs.push_back(s);
    }
    
    sort(seqs.begin(), seqs.end(), cmpSeq);


    for(vector<Seq>::const_iterator it = seqs.begin(); it != seqs.end(); ++it) {
        fs_out << '>' << (*it).desc << endl
               << (*it).data << endl;
    }

    fs_in.close();
    fs_out.close();*/

    return 0;
}
