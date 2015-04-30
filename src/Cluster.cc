#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <thread>

#include "Cluster.h"
#include "Distance.h"
#include "Seq.h"

using namespace std;

/**
 * Given file stream, read sequences of FASTA format, compute the centroids of
 * the clustering (given threshold and k for k-mers) of the given count of
 * sequences. Output the centroids in FASTA format to the output file stream,
 * output clusters to file stream and return the number of centroids.
 */
/*int Cluster::intersect_clust(fstream& fs_in, fstream& fs_centroids, fstream& fs_clusters,
        Distance& dist, int count, int max_rejects) {

    struct Centroid {
        struct Seq seq;
        set<string> kmers;  // set of occurring kmers in seq
    };

    vector<struct Centroid> centroids; // alternatively vector<pair<Seq, set<string>>>

    struct Seq s;
    int i = 0;

    while (IO::read_sequence(fs_in, s) && i < count) {
        bool match = false;
        ++i;

        set<string> query_kmers = dist.kmers(s); // kmers in query sequence

        priority_queue<pair<int, int>> q;  // non-increasingly sorted queue

        // calculate intersection and store cardinality and index in queue
        for (vector<struct Centroid>::const_iterator it = centroids.begin();
                it != centroids.end(); ++it) {
            vector<string> intersect;
            set_intersection(query_kmers.begin(), query_kmers.end(),
                             (it->kmers).begin(), (it->kmers).end(),
                             back_inserter(intersect));
            //cout << "Size of intersect: " << intersect.size() << endl;
            q.push({intersect.size(), it - centroids.begin()});
        }

        int rejects = 0;
        // loop through centroids
        while(!q.empty() && !match && rejects < max_rejects) {
            ++rejects;
            //cout << "q not empty" << endl;
            int ctr_index = (q.top()).second;
            Seq ctr = centroids[ctr_index].seq;
            q.pop();

            if (dist.compare(s.data, ctr.data)) {
                // write s belongs to centroids[i] to fs_clusters
                fs_clusters << ctr_index << ": " << s.data << endl;
                match = true; // found cluster
                break;
            }
        }

        if (!match) {
            // add new centroid and write to stream in FASTA format
            Centroid new_centroid = {s, dist.kmers(s)};
            centroids.push_back(new_centroid);

            fs_centroids << '>' << s.desc << endl
                         << s.data << endl;
        }
    }

    return centroids.size();
}*/


/**
 * Given file stream, read sequences of FASTA format, compute the centroids of
 * the clustering (given threshold and k for k-mers) of the given count of
 * sequences. Output the centroids in FASTA format to the output file stream,
 * output clusters to file stream and return the number of centroids.
 */
int Cluster::simple_clust(const vector<Seq>& seqs, ofstream& fs_centroids,
        ofstream& fs_clusters, Distance& dist, int count, int max_rejects) {
    unordered_map<int, Seq> centroids;

    int centroid_count = 0;

    double write_secs = 0;

    for (auto q_it = seqs.cbegin(); q_it != seqs.cend(); ++q_it) {
        bool match = false;
        vector<int> s_keys = dist.compute_key(*q_it, max_rejects);
        if (s_keys.empty())
          throw logic_error("Calling compute_key on "
                        + (*q_it).to_string() + " returns an empty vector" );

        for (auto it = s_keys.cbegin(); it != s_keys.cend(); ++it) {
            if (centroids.find(*it) == centroids.end()) {
                continue;
            }
            Seq& target = centroids[*it];
            double d = dist.distance(*q_it, target);
            if (d >= dist.threshold()) {
                // write s belongs to centroids[i] to fs_clusters
                //clock_t read_clock = clock();

                /*fs_clusters << "H " << setw(6) << t_it - cts_index.cbegin() << " "
                            << setw(10) << setprecision(5) << fixed << d << " "
                            << (*q_it).desc << " " << seqs[*t_it].desc << "\n";*/

                /*fs_clusters << (*q_it).to_string() << " "
                            << target.to_string()  << '\n';*/

                //write_secs += (clock() - read_clock) / (double) CLOCKS_PER_SEC;

                match = true; // found cluster
                break;
            }
        }

        if (!match) {
            // add new centroid and write to stream in FASTA format
            centroids[s_keys[0]] = *q_it;

            clock_t read_clock = clock();

            // write centroid entry to to clusters file
            fs_clusters << "C " << setw(6) << centroid_count++ << " "
                        << (*q_it).desc << "\n";
            // write FASTA format to centroids file
            fs_centroids << ">" << (*q_it).desc << '\n'
                         << (*q_it).to_string()   << '\n';

            write_secs += (clock() - read_clock) / (double) CLOCKS_PER_SEC;
        }
    }

    cout << "Time spent on file output: " << write_secs << " sec." << endl;

    return centroid_count;
}

int Cluster::thorough_clust(const vector<Seq>& seqs, ofstream& fs_centroids,
        ofstream& fs_clusters, Distance& dist, int count) {

    int centroid_count = 0;
    vector<int> cts_index;

    for (auto q_it = seqs.cbegin(); q_it != seqs.cend(); ++q_it) {
        bool match = false;
        for (auto t_it = cts_index.cbegin(); t_it != cts_index.cend(); ++t_it) {
            double d = dist.distance(*q_it, seqs[*t_it]);
            if (d >= dist.threshold()) {
                fs_clusters << "H " << setw(6) << t_it - cts_index.cbegin() << " "
                            << setw(10) << setprecision(5) << fixed << d << " "
                            << (*q_it).desc << " " << seqs[*t_it].desc << "\n";
                match = true; // found cluster
                break;
            }
        }

        if (!match) {
            // add new centroid and write to stream in FASTA format
            cts_index.push_back(q_it - seqs.cbegin());
            fs_clusters << "C " << setw(6) << centroid_count++ << " "
                        << (*q_it).desc << "\n";
            fs_centroids << (*q_it).to_string() << '\n';
        }
    }
    return centroid_count;
}

/**
 * Fill given bitset with 1's for all the kmers occuring in the given sequence.
 */
inline void get_kmer_bitset(const Seq& s, bitset<KMER_BITSET>& b) {
    static const uint32_t k2 = 2 * KMER_LEN;
    static const uint32_t mask = pow(2, k2) - 1;    // 0b00111..11 (2*k 1's)

    for (size_t i = 0; i <= s.length() - KMER_LEN; ++i) {
        uint32_t kmer = 0;

        //memcpy(&kmer_l, (longer + (i/4)), 4);
        kmer = Distance::stream2int(s.data() + (i/4));
        kmer >>= (32 - 2*(i % 4) - k2);
        kmer &= mask;
        b.set(kmer);
    }
}

/**
 * For every sequence in the given collection, search through the centroids for
 * one where the number of distinct kmers in both the query sequence and
 * centroid are at least dist.threshold()-0.05 times the number of distinct
 * kmers in the centroid sequence. If no match is found after max_rejects
 * tries, the sequence becomes a new centroid.
 */
int Cluster::kmers_select_clust(const vector<Seq>& seqs, ofstream& fs_centroids,
        ofstream& fs_clusters, Distance& dist, int max_rejects) {

    struct Centroid {
        const Seq& seq;
        bitset<KMER_BITSET> bits;
        size_t count;

        Centroid(const Seq& s, bitset<KMER_BITSET> bits) : seq(s), bits(bits) {
            count = bits.count(); // # of distinct kmers in seq
        }
    };

    //typedef bitset<KMER_BITSET> kmer_bits;

    vector<Centroid> centroids;

    unsigned int centroid_count = 0;
    const size_t seqs_size = seqs.size();

    //int neg_count = 0;

    for (auto q_it = seqs.cbegin(); q_it != seqs.cend(); ++q_it) {
        cout << "\r" << 100 * (q_it - seqs.cbegin()) / seqs_size << "%";

        bool match = false;
        int rejects = 0;
        //bool false_negative = false;

        // bitset of kmers occuring in query sequence
        bitset<KMER_BITSET> q_bitset(0);
        get_kmer_bitset(*q_it, q_bitset);

        int i = 0;
        for (auto c_it = centroids.cbegin();
                (c_it != centroids.cend()) && (rejects < max_rejects); ++c_it, ++i) {

            // count number of kmers occurring in both query and target sequence
            size_t set_bits = (q_bitset & c_it->bits).count();

            // if the # of distinct kmers in both query and target is >= to
            // id-0.5 times the # of distinct kmers in the target, then compare
            if (set_bits >= c_it->count * (dist.threshold() - 0.05)) {
                if (dist.compare(*q_it, c_it->seq)) {
                    // write hit entry to to clusters file
                    fs_clusters << 'H' << setw(6) << i
                                << ' ' << (*q_it).desc
                                << ' ' << (c_it->seq).desc;
                    match = true; // found cluster
                    break;
                }
                ++rejects;
            }

            //if (dist.distance(*q_it, c_it->second) >= dist.threshold())
            //    false_negative = true;

        }

        if (!match) {
            //if (false_negative)
            //    ++neg_count;

            // add new centroid and write to stream in FASTA format
            centroids.emplace_back(*q_it, q_bitset);

            /*sort(centroids.begin(), centroids.end(),
                [](const Centroid& a, const Centroid& b) {
                    return a.count > b.count;
                });*/

            // write centroid entry to to clusters file
            fs_clusters << 'C' << setw(6) << centroid_count++ << ' '
                        << (*q_it).desc << '\n';
            // write FASTA format to centroids file
            fs_centroids << '>' << (*q_it).desc << '\n'
                         << (*q_it).to_string() << '\n';

        }
    }

    //cout << endl << "False negatives: " << neg_count << endl;

    cout << "\r100%\n";
    return centroid_count;
}

void sub_clust(vector<Seq>::iterator begin, vector<Seq>::iterator end,
        vector<Centroid>& cts, Distance& dist, int max_rejects) {
    //const size_t seqs_size = end - begin;
    const size_t seqs_size = distance(begin, end);

    for (auto q_it = begin; q_it != end; ++q_it) {
        cout << "\r" << 100 * (q_it - begin) / seqs_size << "%";

        bool match = false;
        int rejects = 0;

        // bitset of kmers occuring in query sequence
        bitset<KMER_BITSET> q_bitset(0);
        get_kmer_bitset(*q_it, q_bitset);

        int i = 0;
        for (auto c_it = cts.begin();
                (c_it != cts.end()) && (rejects < max_rejects); ++c_it, ++i) {
            // count number of kmers occurring in both query and target sequence
            size_t set_bits = (q_bitset & c_it->bits).count();

            // if the # of distinct kmers in both query and target is >= to
            // id-0.5 times the # of distinct kmers in the target, then compare
            if (set_bits >= c_it->count * (dist.threshold() - 0.05)) {
                if (dist.compare(*q_it, c_it->seq)) {
                    (c_it->cls_seqs).push_back(ref(*q_it));
                    //(c_it->cls_seqs).push_back(move(*q_it)); // TODO: maybe move?
                    match = true; // found cluster
                    break;
                }
                ++rejects;
            }
        }

        if (!match) {
            // add new centroid and write to stream in FASTA format
            cts.emplace_back(*q_it, move(q_bitset)); // TODO: does the move work as expected?
        }
    }
    cout << "\r100%\n";
    //return cts.size();
}


int Cluster::clust(vector<Seq>& seqs, Distance& dist, int max_rejects) {

    cout << thread::hardware_concurrency()
         << " concurrent threads are supported.\n"; // only a hint

    unsigned int mid = (seqs.end() - seqs.begin()) / 2;

    vector<Seq>::iterator fst_beg = seqs.begin();
    vector<Seq>::iterator fst_end = seqs.begin() + mid;
    vector<Seq>::iterator snd_beg = seqs.begin() + mid + 1;
    vector<Seq>::iterator snd_end = seqs.end();

    vector<Centroid> c0, c1;

    thread t0(sub_clust, fst_beg, fst_end, ref(c0), ref(dist), max_rejects);
    thread t1(sub_clust, snd_beg, snd_end, ref(c1), ref(dist), max_rejects);

    t0.join();
    cout << "t0 joined" << endl;
    t1.join();
    cout << "t1 joined" << endl;

    /*sub_clust(seqs.begin(), seqs.begin() + mid, c0, dist, max_rejects);
    sub_clust(seqs.begin() + mid + 1, seqs.end(), c1, dist, max_rejects);*/

    cout << "Finished divide. Number of clusters: "
         << c0.size() << ", " << c1.size() << "\n"
         << "Combining clusters...\n";

    int rejects = 0;

    auto c0_end = c0.end();

    for (auto it1 = c1.begin(); it1 != c1.end(); ++it1) {
        bool match = false;

        for (auto it0 = c0.begin();
                it0 != c0_end && rejects < max_rejects; ++it0) {

            size_t set_bits = (it0->bits & it1->bits).count();
            if (set_bits >= it0->count * (dist.threshold() - 0.05)) {
                if (dist.compare(it0->seq, it1->seq)) {
                    // merge clusters, i.e. combine vectors of sequences
                    (it0->cls_seqs).insert((it0->cls_seqs).end(),
                            (it1->cls_seqs).begin(), (it1->cls_seqs).end());
                    //(it0->cls_seqs).push_back(it1->seq);
                    match = true;
                    break;
                }
                ++rejects;
            }
        }

        if (!match) {
            c0.push_back(*it1);
        }
    }

    //vector<Centroid>::iterator cts0_end = cts0.end();
    /*for (auto c1 : cts1) {
        for (auto c0 : cts0) {
            if (dist.compare(c0.seq, c1.seq))
                break;
            ++centroids_count;
        }
    }*/

    return c0.size();
}
