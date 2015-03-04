#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "Distance.h"

//using namespace std;

Distance dist;

/**
 * Read [DR]NA sequence from given filestream in given string.
 * Return 0 on success and -1 if there's no more to read.
 */
int readSequence(std::fstream& fs, std::string& s) {
    std::string tmp;
    s = ""; // clear s

    std::getline(fs, tmp);
    if (tmp.empty())
        return -1;
    while (tmp[0] == '>')
        std::getline(fs, tmp);                  // ignore '>' line
    while (!tmp.empty() && tmp[0] != '>') {     // read until next '>' line
        s += tmp;
        std::getline(fs, tmp);
    }
    return 0;
}


int main(int argc, char *argv[])
{
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <name of .fasta file> "
                                             " <k in k-mers> "
                                             " <# of sequences to compare>"
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

    int distances[count][count];

    for (int i = 0; i < count; i++) {
        for (int j = 0; j < count; j++) {
            distances[i][j] = -1;
        }
    }

    for (int i = 0; i < count; ++i) {
        readSequence(fs0, fst);

        for (int j = 0; j < count; ++j) {
            readSequence(fs1, snd);
            if (i == j) { // don't compare a sequence to itself
                distances[i][j] = 0;
                continue;
            }
            if (distances[i][j] != -1) {
                continue;
            }
            int newdist = dist.d2window(fst, snd, k);
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

/*    std::string a = "ctg";
    std::string b = "caacat";
    std::cout << dist.d2window(a,b,2);
    return 0;*/
}