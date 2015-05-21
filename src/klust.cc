#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <list>
#include <string>
#include <unistd.h>

#include "Cluster.h"
#include "Distance.h"
#include "IO.h"
#include "Seq.h"
#include "Utils.h"

using namespace std;

int main(int argc, char *argv[])
{

    ios_base::sync_with_stdio(false); // don't share buffers with C style IO

    // default similarity and clustering parameters
    int k = 6;
    double thrs = 0.85;
    int count = INT_MAX;
    int max_rejects = 8;
    int step = 1;
    bool sort_incr = false;
    bool sort_decr = false;
    int depth = 0;

    // misc parameters
    bool springy = false;

    // CLI argument parsing
    static struct option long_options[] = {
        {"count",       required_argument, 0, 'c'},
        {"depth",       required_argument, 0, 'l'},
        {"id",          required_argument, 0, 't'},
        {"max_rejects", required_argument, 0, 'm'},
        {"sort_decr",   no_argument,       0, 'd'},
        {"sort_incr",   no_argument,       0, 'i'},
        {"step_size",   required_argument, 0, 's'},

        {"springy",     no_argument,       0,  0 }
    };

    int opt, option_index = 0;
    while ((opt = getopt_long(argc, argv, "c:k:l:m:s:t:",
                    long_options, &option_index)) != -1) {
        switch (opt) {
            case 'c':
                count = atoi(optarg);
                break;
            case 'd':
                sort_decr = true;
                break;
            case 'i':
                sort_incr = true;
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'l':
                depth = atoi(optarg);
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

            // misc other parameters
            case  0 :
                springy = true;
                break;
            default:
                cout << "unexpected argument" << endl;
                return 1;
        }
    }

    if (argc < (optind + 3)) {
        cerr << "Usage: " << argv[0] << " <FASTA input file> "
                                        " <FASTA output file for centroid> "
                                        " <output file for clustering results> "
                                     << endl; // TODO: mention options
        return 1;
    }

    if (sort_incr && sort_decr) {
        cerr << "Error: can't sort both by both "
                "decreasing and increasing length." << endl;
        return 1;
    }

    ifstream fs_in(argv[optind++]);
    ofstream fs_cts(argv[optind++], ofstream::out | ofstream::trunc);
    ofstream fs_cls(argv[optind++], ofstream::out | ofstream::trunc);

    if (!(fs_in && fs_cts && fs_cls)) {
        cerr << "Error opening file" << endl;
        fs_in.close();
        fs_cts.close();
        fs_cls.close();
        return 1;
    }

    cout << "Running with parameters: "
         << "\n  k = " << k
         << "\n  id = " << thrs
         << "\n  max_rejects = " << max_rejects
         << "\n  depth = " << depth
         << "\n" << endl;

    /*
     * Reading sequences
     */
    vector<Seq> seqs;
    cout << "Reading sequences..." << endl;

    // timing
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    count = IO::read_seqs(fs_in, seqs, count);

    end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;

    cout << "Time: "     << elapsed_seconds.count() << " sec.\n"
         << "Seqs/sec: " << count / elapsed_seconds.count() << "\n" << endl;

    if (sort_decr) {
        cout << "Sorting by decreasing sequence length..." << endl;
        sort(seqs.begin(), seqs.end(),
            [](Seq& s1, Seq& s2) {
                return s1.length > s2.length;
            });
    }
    if (sort_incr) {
        cout << "Sorting by increasing sequence length..." << endl;
        sort(seqs.begin(), seqs.end(),
            [](Seq& s1, Seq& s2) {
                return s1.length < s2.length;
            });
    }

    /*Utils::permute(seqs, 10, 0.01, fs_cts);
    return 0;*/
    /*Utils::permute_chunks(seqs, 10, 0.01, fs_cts, 5);
    return 0;*/

 /*   Distance dist(k, thrs, step);
    Utils::print_matrix(seqs, cout, dist);
    return 0;*/
    


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
    Cluster clust(Distance(k, thrs, step), max_rejects);
    list<Centroid> cts;

    cout << "Clustering " << count << " sequences..." << endl;

    // timing
    start = chrono::system_clock::now();

    clust.clust(seqs, cts, depth);

    end = chrono::system_clock::now();
    elapsed_seconds = end - start;

    cout << '\n'
         << "Time: "        << elapsed_seconds.count() << " sec.\n"
         << "Throughput: "  << count / elapsed_seconds.count() << " seqs/sec.\n"
         << endl;

    /*
     * Outputting results
     */
    IO::print_stats(seqs, cts);

    IO::write_results(cts, fs_cts, fs_cls);

    if (springy) {
        ofstream fs_springy("springy.html", ofstream::out | ofstream::trunc);
        if (fs_springy) {
            IO::springy(cts, fs_springy);
            fs_springy.close();
        } else
            cerr << "error generating springy code" << endl;
    }

    fs_in.close();
    fs_cts.close();
    fs_cls.close();

    return 0;
}
