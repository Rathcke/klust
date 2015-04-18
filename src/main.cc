#include <algorithm>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <unistd.h>

#include "Cluster.h"
#include "Distance.h"
#include "IO.h"

using namespace std;

int main(int argc, char *argv[])
{
    ios_base::sync_with_stdio(false); // don't share buffers with C style IO

    // default similarity and clustering parameters
    int k = 6;
    double thrs = 0.85;
    int count = INT_MAX; // SILVA: 1583830, TODO: INT_MAX maybe not so pretty.
    int max_rejects = 8;
    int step = 1;

    // CLI argument parsing
    static struct option long_options[] = {
        {"count",       required_argument, 0, 'c'},
        {"id",          required_argument, 0, 't'},
        {"max_rejects", required_argument, 0, 'm'},
        {"step_size",   required_argument, 0, 's'}
    };

    int opt, option_index = 0;
    while ((opt = getopt_long(argc, argv, "c:k:m:s:t:",
                    long_options, &option_index)) != -1) {
        switch (opt) {
            case 'c':
                count = atoi(optarg);
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'm':
                max_rejects = atoi(optarg);
                break;
            case 's':
                step = atoi(optarg);
                break;
            case 't':
                thrs = atof(optarg);
                break;
            default:
                cout << "unexpected argument" << endl;
                return 1;
        }
    }

    if ((optind + 3) > argc) {
        cerr << "Usage: " << argv[0] << " <FASTA input file> "
                                        " <FASTA output file for centroid> "
                                        " <output file for clustering results> "
                                     << endl; // TODO: mention options
        return 1;
    }

    ifstream fs_in(argv[optind++]);
    ofstream fs_cts(argv[optind++], ofstream::out | ofstream::trunc);
    ofstream fs_cls(argv[optind++], ofstream::out | ofstream::trunc);

    if (!(fs_in && fs_cts && fs_cls)) {
        cerr << "error opening file" << endl;
        fs_in.close();
        fs_cts.close();
        fs_cls.close();
        return 1;
    }

    cout << "Runnig with parameters: "
         << "\n  k = " << k
         << "\n  id = " << thrs
         << "\n  count = " << count
         << "\n  max_rejects = " << max_rejects
         << "\n  step_size = " << step
         << "\n" << endl;
         
	Distance d2(k, thrs, step);

    /*
     * Reading sequences
     */
    vector<Seq> seqs;

    cout << "Reading " << count << " sequences..." << endl;
    clock_t read_clock = clock();
    IO::read_seqs(fs_in, seqs, count);
    double read_secs = (clock() - read_clock) / (double) CLOCKS_PER_SEC;

    cout << "Finished reading:\n"
         << "Time: "     << read_secs << " sec.\n"
         << "Seqs/sec: " << count / read_secs << "\n" << endl;

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
    cout << "# of clusters: "
         << Cluster::simple_clust(seqs, fs_cts, fs_cls, d2, count, max_rejects)
         //<< Cluster::thorough_clust(seqs, fs_cts, fs_cls, d2, count)
         << endl;
    double comp_secs = (clock() - comp_clock) / (double) CLOCKS_PER_SEC;

    cout << "Finished clustering:\n"
         << "Time: "            << comp_secs << " sec.\n"
         << "Throughput: "      << count / comp_secs << " seqs/sec.\n";

    fs_in.close();
    fs_cts.close();
    fs_cls.close();

    return 0;
}
