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
    ios_base::sync_with_stdio(false); // don't share buffers with C style IO

    ifstream fs_in(argv[1]);
    ofstream fs_cts(argv[2], ofstream::out | ofstream::trunc);
    ofstream fs_cls(argv[3], ofstream::out | ofstream::trunc);

    if (!fs_in || !fs_cts || !fs_cls) {
        cerr << "error opening file: " << argv[1] << endl;
        return 1;
    }

	Distance d2(6, 0.85, 0);
    int count = 100 * 1000;
    int max_rejects = 8;

    /*
     * Reading sequences
     */
    /*int count = 500;
    //vector<vector<bitset<2>>> seqs;
    vector<Seq> seqs;

    cout << "Reading " << count << " sequences...\n" << endl;
    clock_t read_clock = clock();
    IO::read_seqs(fs_in, seqs, count);
    double read_secs = (clock() - read_clock) / (double) CLOCKS_PER_SEC;

    cout << "Finished reading:\n"
         << "Time: "     << read_secs << " sec.\n"
         << "Seqs/sec: " << count / read_secs << "\n" << endl;*/

    /*
     * Comparing sequences
     */
    /*int tot = 0;
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
    cout << tot << endl;*/

    /*
     * Clustering
     */
    cout << "Clustering " << count << " sequences..." << endl;
    clock_t comp_clock = clock();
    cout << "# of clusters: " <<
        Cluster::simple_clust(fs_in, fs_cts, fs_cls, d2, count, max_rejects) << endl;
    double comp_secs = (clock() - comp_clock) / (double) CLOCKS_PER_SEC;

    cout << "Finished clustering:\n"
         << "Time: "            << comp_secs << " sec.\n"
         << "Throughput: "      << count / comp_secs << "seqs/sec.\n";

    fs_in.close();
    fs_cts.close();
    fs_cls.close();

    return 0;
}
