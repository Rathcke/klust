#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <assert.h>

#include "IO.h"
#include "Distance.h"
#include "Cluster.h"

//using namespace std;

IO io;
Distance dist;
Cluster cluster;

int main(int argc, char *argv[])
{
    if (argc < 5) {
        std::cout << "Usage: " << argv[0] << " <name of .fasta file> "
                                             " <k in k-mers> "
                                             " <# of sequences to compare>"
                                             " <Similarity threshold>"
                                          << std::endl << std::endl;
        std::cout << "This program will measure the d2 distance between \n"
                     "all of the specified number of sequences from the \n"
                     "given .fasta file and output a distance matrix.\n"
                  << std::endl;
        return 1;
    }

    std::fstream fs0(argv[1]);
    std::fstream fs1(argv[1]);

    std::string fst, snd;

    const int k = std::atoi(argv[2]);       // k in k-mer
    const int count = std::atoi(argv[3]);   // # of sequences to measure
    const int threshold = std::atoi(argv[4]); // Simlilarity threshold

    std::cout << cluster.clust(fs0, fs1, threshold, k, count) << std::endl;

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

    return 0;
}
