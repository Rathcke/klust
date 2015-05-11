#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

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

    for (size_t i = 0; i <= s.length - KMER_LEN; ++i) {
        uint32_t kmer = 0;

        //memcpy(&kmer_l, (longer + (i/4)), 4);
        kmer = Distance::stream2int(s.data + (i/4));
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
void Cluster::kmer_select_clust(vector<Seq>::const_iterator begin,
        vector<Seq>::const_iterator end, list<Centroid>& cts) {
    const size_t seqs_size = distance(begin, end);
    unsigned int centroid_count = 0;

    for (auto q_it = begin; q_it != end; ++q_it) {
        cout << "\r" << 100 * (q_it - begin) / seqs_size << "%";

        bool match = false;
        int rejects = 0;

        // bitset of kmers occuring in query sequence
        bitset<KMER_BITSET> q_bitset(0);
        get_kmer_bitset(*q_it, q_bitset);

        const Seq *close_match = nullptr;

        int i = 0;
        for (auto c_it = cts.begin();
                (c_it != cts.end()) && (rejects < max_rejects); ++c_it, ++i) {
            // count number of kmers occurring in both query and target sequence
            size_t set_bits = (q_bitset & c_it->bits).count();

            // if the # of distinct kmers in both query and target is >= to
            // id-0.5 times the # of distinct kmers in the target, then compare
            if (set_bits >= c_it->count * dist.threshold()) {
                if (dist.compare(*q_it, c_it->seq)) {
                    (c_it->cls_seqs).push_back(ref(*q_it));
                    match = true; // found cluster
                    cts.push_front(move(*c_it));
                    cts.erase(c_it);
                    break;
                }
                if (c_it->link) {
                    if (dist.compare(*q_it, *(c_it->link))) {
                        (c_it->cls_seqs).push_back(ref(*q_it));
                        match = true; // found cluster
                        break;
                    }
                }
                close_match = &(c_it->seq);
                ++rejects;
            }
        }

        if (!match) {
            // add new centroid to list
            cts.emplace_front(*q_it, q_bitset, centroid_count++);
            if (close_match)
                cts.front().link = close_match;
        }
    }
    cout << "\r100%";
}


int Cluster::clust(const vector<Seq>& seqs, list<Centroid>& cts, int depth) {
    cout << thread::hardware_concurrency()
         << " concurrent threads are supported.\n"; // only a hint

    int seqs_size = seqs.size();
    int sub_count = pow(2, depth);
    int sub_size = seqs_size / sub_count;

    vector<list<Centroid>> cts_ls;
    cts_ls.resize(sub_count);

    vector<thread> threads;

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
            ref(cts_ls[sub_count-1]));

    for (auto& t : threads)
        t.join();

    int clusters_sum = 0;
    cout << "\nFinished divide. Number of clusters: ";
    for (auto& ls : cts_ls) {
        int ls_size = ls.size();
        clusters_sum += ls_size;
        cout << ls_size << " ";
    }
    cout << "\nCombining clusters..." << endl;

    // merge centroid lists in a bottom-up manner; result will be in first list
    int delta = 2;
    int i = 0;
    while (i < sub_count-1 && delta <= sub_count) {
        //cout << "merging centroids " << i << " and " << i+(delta/2) << endl;
        merge(cts_ls[i], cts_ls[i + (delta/2)]);
        if (i + delta >= sub_count-1) {
            i = 0;
            delta *= 2;
        } else {
            i += delta;
        }
    }

    cout << "Merged " << clusters_sum - cts_ls[0].size() << " clusters." << endl;

    cts = move(cts_ls[0]); // move resulting centroid list into argument
    return cts.size();
}

/**
 * Merge two vectors of centroids and store the result in the first vector.
 */
void Cluster::merge(list<Centroid>& res, const list<Centroid>& c1) {
    int rejects = 0;
    const auto res_end = res.end();

    for (auto it1 = c1.begin(); it1 != c1.end(); ++it1) {
        bool match = false;

        for (auto it0 = res.begin();
                it0 != res_end && rejects < max_rejects; ++it0) {

            size_t set_bits = (it0->bits & it1->bits).count();
            if (set_bits >= it0->count * (dist.threshold() - 0.05)) {
                if (dist.compare(it0->seq, it1->seq)) {
                    // merge clusters, i.e. combine vectors of sequences
                    (it0->cls_seqs).insert((it0->cls_seqs).end(),
                            (it1->cls_seqs).begin(), (it1->cls_seqs).end());
                    (it0->cls_seqs).push_back(it1->seq); // push_back centroid
                    match = true;
                    break;
                }
                ++rejects;
            }
        }
        if (!match)
            res.push_back(*it1);
    }
}
