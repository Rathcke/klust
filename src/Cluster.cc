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

int Cluster::simple_clust(const vector<Seq>& seqs, ofstream& fs_centroids,
        ofstream& fs_clusters) {

    unordered_map<int, Seq> centroids;

    int centroid_count = 0;

    for (auto q_it = seqs.cbegin(); q_it != seqs.cend(); ++q_it) {

        bool match = false;
        vector<int> s_keys = dist.compute_key(*q_it, max_rejects);

        if (s_keys.empty())
          throw logic_error("Calling compute_key on "
                  + (*q_it).to_string() + " returns an empty vector");

        for (auto it = s_keys.cbegin(); it != s_keys.cend(); ++it) {
            if (centroids.find(*it) == centroids.end())
                continue;

            Seq& target = centroids[*it];

            if (dist.compare(*q_it, target)) {
                // write s belongs to centroids[i] to fs_clusters
                /*fs_clusters << "H " << setw(6) << t_it - cts_index.cbegin() << " "
                            << setw(10) << setprecision(5) << fixed << d << " "
                            << (*q_it).desc << " " << seqs[*t_it].desc << "\n";

                fs_clusters << (*q_it).to_string() << " "
                            << target.to_string()  << '\n';*/

                match = true; // found cluster
                break;
            }
        }

        if (!match) {
            // add new centroid and write to stream in FASTA format
            centroids[s_keys[0]] = *q_it;

            // write centroid entry to to clusters file
            fs_clusters << "C " << setw(6) << centroid_count++ << " "
                        << (*q_it).desc << "\n";
            // write FASTA format to centroids file
            fs_centroids << ">" << (*q_it).desc << '\n'
                         << (*q_it).to_string()   << '\n';
        }
    }

    return centroid_count;
}

int Cluster::thorough_clust(const vector<Seq>& seqs, ofstream& fs_centroids,
        ofstream& fs_clusters) {

    int centroid_count = 0;
    vector<int> cts_index;
    int numb = -1;

    for (auto q_it = seqs.cbegin(); q_it != seqs.cend(); ++q_it) {
        ++numb;
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
            cout << numb << ", ";
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

        kmer = Distance::stream2int(s.data + (i/4));
        kmer >>= (32 - 2*(i % 4) - k2);
        kmer &= mask;
        b.set(kmer);
    }
}

void Cluster::kmer_select_clust(vector<Seq>::const_iterator begin,
        vector<Seq>::const_iterator end, list<Centroid>& cts) {

    const size_t seqs_size = distance(begin, end);
    unsigned int centroid_count = 0;

    bitset<KMER_BITSET> q_bitset(0);

    for (auto q_it = begin; q_it != end; ++q_it) {
        cout << "\r" << 100 * (q_it - begin) / seqs_size << "%";

        bool match = false;
        int rejects = 0;    // number of unsuccessful compares so far

        // bitset of kmers occuring in query sequence
        q_bitset.reset();
        get_kmer_bitset(*q_it, q_bitset);

        // pointer to most recent unsuccesful, compared centroid Seq
        Centroid *close_match = nullptr;

        for (auto c_it = cts.begin();
                (c_it != cts.end()) && (rejects < max_rejects); ++c_it) {

            // count number of kmers occurring in both query and target sequence
            size_t set_bits = (q_bitset & c_it->bits).count();

            // if the # of distinct kmers in both query and target is >= to
            // id times the # of distinct kmers in the target, then compare
            if (set_bits >= c_it->count * dist.threshold()) {
                if (dist.compare(*q_it, c_it->seq)) {
                    (c_it->cls_seqs).push_back(ref(*q_it));
                    match = true; // found cluster
                    // moves element at c_it to front of cts
                    cts.splice(cts.begin(), cts, c_it);
                    break;
                }
                if (c_it->link) {
                    if (dist.compare(*q_it, c_it->link->seq)) {
                        (c_it->link->cls_seqs).push_back(ref(*q_it));
                        match = true; // found cluster
                        break;
                    }
                }
                ++rejects;
                close_match = &(*c_it);
            }
        }

        if (!match) {
            // add new centroid to list
            cts.emplace_front(*q_it, q_bitset, centroid_count++);
            //cout << numb << ", ";
            if (close_match)
                cts.front().link = close_match;
        }
    }

    cout << "\r100%";
}


int Cluster::clust(const vector<Seq>& seqs, list<Centroid>& cts, int depth) {
    if (depth == 0) {
        kmer_select_clust(seqs.begin(), seqs.end(), cts);
        return cts.size();
    }

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

void Cluster::merge(list<Centroid>& res, const list<Centroid>& c1) {
    int rejects = 0;
    const auto res_end = res.end();

    for (auto it1 = c1.begin(); it1 != c1.end(); ++it1) {
        bool match = false;

        for (auto it0 = res.begin();
                it0 != res_end && rejects < max_rejects; ++it0) {

            size_t set_bits = (it0->bits & it1->bits).count();
            if (set_bits >= it0->count * dist.threshold()) {
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
