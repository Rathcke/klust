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

IO io;
Distance dist;
Cluster cluster;

int main(int argc, char *argv[])
{
    if (argc < 6) {
        std::cout << "Usage: " << argv[0] << " <.fasta input file> "
                                             " <.fasta output file for centroids> "
                                             " <k in k-mers> "
                                             " <# of sequences to compare>"
                                             " <similarity threshold>"
                                          << std::endl << std::endl;
        std::cout << std::endl;
        return 1;
    }

    std::fstream fs0(argv[1], fstream::in);
    std::fstream fs1(argv[1], fstream::in);
    std::fstream fs2(argv[2],  // new file doesn't open witout trunc flag
            fstream::in | fstream::out | fstream::trunc);

    const int k = std::atoi(argv[3]);         // k in k-mer
    const int count = std::atoi(argv[4]);     // # of sequences to measure
    const int threshold = std::atoi(argv[5]); // Simlilarity threshold

    cout << "# of clusters: " <<
        cluster.clust(fs0, fs2, threshold, k, count) << endl;


    std::string fst, snd;

    int distances[count][count];
    
    fs0.seekg(0, std::ios::beg);

    for (int i = 0; i < count; i++) {
        for (int j = 0; j < count; j++) {
            distances[i][j] = -1;
        }
    }

    for (int i = 0; i < count; ++i) {
        io.readSequence(fs0, fst);

        for (int j = 0; j < count; ++j) {
            io.readSequence(fs1, snd);
            if (i == j) { // don't compare a sequence to itself
                distances[i][j] = 0;
                continue;
            }
            if (distances[i][j] != -1) {
                continue;
            }

            int newdist = dist.d2window(fst, snd, k);

            // TODO: testing equal to naive version
            //assert(newdist == dist.d2window_naive(fst, snd, k));
            distances[i][j] = newdist;
            distances[j][i] = newdist;
        }
        fs1.seekg(0, std::ios::beg); // rewind fs1 to start
    }

    for (int i = 0; i < count; ++i) {
        //std::cout << "i = " << std::setw(3) <<  i << " : ";
        for (int j = 0; j < count; ++j) {
            std::cout << std::setw(4) << distances[i][j];
        }
        std::cout << std::endl;
    }

    fs0.close();
    fs1.close();
    fs2.close();

    return 0;
}
