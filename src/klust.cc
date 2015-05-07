#include <algorithm>
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

#include <thread>

#include "Cluster.h"
#include "Distance.h"
#include "IO.h"

using namespace std;

/**
 * Fill given bitset with 1's for all the kmers occuring in the given sequence.
 */
void set_kmer_bitset(vector<Seq>::iterator b, vector<Seq>::iterator e) {
    static const uint32_t k2 = 2 * KMER_LEN;
    static const uint32_t mask = pow(2, k2) - 1;    // 0b00111..11 (2*k 1's)

    bitset<KMER_BITSET> bits;

    for (auto it = b; it != e; ++it) {
        for (size_t i = 0; i <= (*it).length() - KMER_LEN; ++i) {
            uint32_t kmer = 0;

            //memcpy(&kmer_l, (longer + (i/4)), 4);
            kmer = Distance::stream2int((*it).data() + (i/4));
            kmer >>= (32 - 2*(i % 4) - k2);
            kmer &= mask;
            bits.set(kmer);
        }
        (*it).set_bits(bits);
    }
}

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

    // CLI argument parsing
    static struct option long_options[] = {
        {"count",       required_argument, 0, 'c'},
        {"depth",       required_argument, 0, 'l'},
        {"id",          required_argument, 0, 't'},
        {"max_rejects", required_argument, 0, 'm'},
        {"sort_decr",   no_argument,       0, 'd'},
        {"sort_incr",   no_argument,       0, 'i'},
        {"step_size",   required_argument, 0, 's'}
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
    clock_t read_clock = clock();
    count = IO::read_seqs(fs_in, seqs, count);
    double read_secs = (clock() - read_clock) / (double) CLOCKS_PER_SEC;
    cout << "Time: "     << read_secs << " sec.\n"
         << "Seqs/sec: " << count / read_secs << "\n" << endl;

    if (sort_decr) {
        cout << "Sorting by decreasing sequence length..." << endl;
        sort(seqs.begin(), seqs.end(),
            [](Seq& s1, Seq& s2) {
                return s1.length() > s2.length();
            });
    }
    if (sort_incr) {
        cout << "Sorting by increasing sequence length..." << endl;
        sort(seqs.begin(), seqs.end(),
            [](Seq& s1, Seq& s2) {
                return s1.length() < s2.length();
            });
    }

    vector<thread> threads;

    int sub_count = 1; // thread::hardware_concurrency();
    int seqs_size = seqs.size();
    int sub_size = seqs_size / sub_count;

    cout << "Calculating bitsets using " << sub_count << " threads..." << endl;
    for (int i = 0; i < sub_count; ++i) {
        threads.emplace_back(set_kmer_bitset, seqs.begin() + i*sub_size,
                seqs.begin() + (i+1)*sub_size);
        cout << (seqs.begin() + i*sub_size) - seqs.begin() << " : ";
        cout << (seqs.begin() + (i+1)*sub_size) - seqs.begin() << endl;
    }
    for (auto& t : threads)
        t.join();
    cout << "Finished!" << endl;

    //cout << seqs[0].get_bits() << endl;
    for (auto& s : seqs)
        cout << s.get_count() << " ";
    cout << endl;
    //cout << seqs[0].get_count() << endl;


    /*vector<list<Centroid>> cts_ls;
    cts_ls.resize(sub_count);


    cout << "Diving " << depth << " times ("
         << sub_count << " sub-clustering(s) and "
         << sub_count << " concurrent thread(s))..." << endl;
    for (int i = 0; i < sub_count-1; ++i) {
        threads.emplace_back(&Cluster::kmer_select_clust, this,
                seqs.begin() + i*sub_size, seqs.begin() + (i+1)*sub_size - 1,
                ref(cts_ls[i]));
    }
    threads.emplace_back(&Cluster::kmer_select_clust, this,
            seqs.begin() + sub_size*(sub_count-1), seqs.end(),
            ref(cts_ls[sub_count-1]));*/


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

    cout << "Kmers Select Clustering " << count << " sequences..." << endl;
    clock_t comp_clock = clock();
    cout << "# of clusters: "
         << clust.clust(seqs, cts, depth)
         << endl;
    double comp_secs = (clock() - comp_clock) / (double) CLOCKS_PER_SEC;

    cout << "Finished clustering:\n"
         << "Time: "            << comp_secs << " sec.\n"
         << "Throughput: "      << count / comp_secs << " seqs/sec.\n";

    /*
     * Outputting results
     */
    IO::write_results(cts, fs_cts, fs_cls);


    fs_in.close();
    fs_cts.close();
    fs_cls.close();

    return 0;
}
