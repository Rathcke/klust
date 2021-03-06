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

void print_usage(char *argv[]);

int main(int argc, char *argv[])
{
    ios_base::sync_with_stdio(false); // don't share buffers with C style IO

    // default similarity and clustering parameters
    int k = 5;
    double thrs = 0.85;
    int count = INT_MAX;
    int max_rejects = 8;
    bool sort_incr = false;
    bool sort_decr = false;
    int depth = 0;

    bool output_centroids = false;
    bool output_clusters = false;
    string centroids_file, clusters_file;

    // misc parameters
    bool springy = false;
    string springy_file;

    // CLI argument parsing
    static struct option long_options[] = {
        {"centroids",   required_argument, 0, 'o'},
        {"clusters",    required_argument, 0, 'u'},
        {"count",       required_argument, 0, 'c'},
        {"depth",       required_argument, 0, 'l'},
        {"id",          required_argument, 0, 't'},
        {"max_rejects", required_argument, 0, 'm'},
        {"sort_decr",   no_argument,       0, 'd'},
        {"sort_incr",   no_argument,       0, 'i'},

        {"springy",     required_argument, 0,  1 },
        {0,             0,                 0,  0 }
    };

    int opt, option_index = 0;
    while ((opt = getopt_long(argc, argv, "c:dik:l:m:o:t:u:",
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
            case 'o':
                output_centroids = true;
                centroids_file = optarg;
                break;
            case 't':
                thrs = atof(optarg);
                break;
            case 'u':
                output_clusters = true;
                clusters_file = optarg;
                break;

            // misc other parameters
            case  1 :
                springy = true;
                springy_file = optarg;
                break;
            default:
                cout << "unexpected argument" << endl;
                print_usage(argv);
                return 1;
        }
    }

    if (argc < optind + 1) {
        print_usage(argv);
        return 0;
    }

    if (sort_incr && sort_decr) {
        cerr << "Error: can't sort both by both "
                "decreasing and increasing length." << endl;
        return 1;
    }

    ifstream fs_in(argv[optind]);
    if (!fs_in) {
        cerr << "Error opening file \"" << argv[optind] << "\"" << endl;
        return 1;
    }
    ++optind;

    ofstream fs_cts, fs_cls;

    if (output_centroids) {
        fs_cts.open(centroids_file, ofstream::out | ofstream::trunc);
        if (!fs_cts) {
            cerr << "Error opening file \"" << centroids_file << "\"" << endl;
            fs_in.close();
            return 1;
        }
    }

    if (output_clusters) {
        fs_cls.open(clusters_file, ofstream::out | ofstream::trunc);
        if (!fs_cls) {
            cerr << "Error opening file \"" << clusters_file << "\"" << endl;
            fs_in.close();
            fs_cts.close();
            return 1;
        }
    }

    cout << "Running with parameters: "
         << "\n  k = " << k
         << "\n  id = " << thrs
         << "\n  max_rejects = " << max_rejects
         << "\n  depth = " << depth
         << "\n" << endl;

    /*
     * Setting up distance metric
     */
    Distance dist(k, thrs);

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

    /*
     * Clustering
     */
    Cluster clust(dist, max_rejects);
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

    if (output_centroids)
        IO::write_centroids(cts, fs_cts);
    if (output_clusters)
        IO::write_clusters(cts, fs_cls);

    if (springy) {
        ofstream fs_springy(springy_file, ofstream::out | ofstream::trunc);
        if (fs_springy) {
            cout << "Generating springy code..." << endl;
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

void print_usage(char *argv[]) {
    cout << "Usage: " << argv[0] << " <FASTA input file>\n"
            "\n"
            "Options:\n"
            "-o, --centroids file   Output FASTA file for centroids\n"
            "-u, --clusters file    Output file for clustering results\n"
            "-t, --id t             Set similarity threshold/identity to t (in [0,1])\n"
            "-k k                   Set the k in k-mer, used in similarity metric\n"
            "-c, --count n          Read and process n sequences\n"
            "-d, --sort_decr        Sort sequences by decresing length\n"
            "-i, --sort_incr        Sort sequences by increasing length\n"
            "-l, --depth d          Set the depth of the tree of divides, i.e.\n"
            "                       cluster 2^d subsets of sequences and combine\n"
            "-m, --max_rejects x    Max number of rejects when searching for centroid\n"
            "\n"
            "--springy file         Generate springy JavaScript code\n"
            "\n"
            "Example: \n"
            "./klust ../data/SILVA_10k.fasta -o centroids.fasta -u clusters --sort_incr\n";
}
