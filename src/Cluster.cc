#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "Cluster.h"
#include "Distance.h"
#include "IO.h"

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
                            << (*q_it).desc() << " " << seqs[*t_it].desc() << "\n";*/

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
                        << (*q_it).desc() << "\n";
            // write FASTA format to centroids file
            fs_centroids << ">" << (*q_it).desc() << '\n'
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
                            << (*q_it).desc() << " " << seqs[*t_it].desc() << "\n";
                match = true; // found cluster
                break;
            }
        }

        if (!match) {
            // add new centroid and write to stream in FASTA format
            cts_index.push_back(q_it - seqs.cbegin());
            fs_clusters << "C " << setw(6) << centroid_count++ << " "
                        << (*q_it).desc() << "\n";
            fs_centroids << (*q_it).to_string() << '\n';
        }
    }
    return centroid_count;
}

int Cluster::kmers_select_clust(const vector<Seq>& seqs, ofstream& fs_centroids,
        ofstream& fs_clusters, Distance& dist, int max_rejects) {

    typedef bitset<KMER_BITSET> kmer_bits;
    vector<pair<kmer_bits, Seq>> centroids;

    int centroid_count = 0;
    size_t count = seqs.cend() - seqs.cbegin();

    //int neg_count = 0;

    for (auto q_it = seqs.cbegin(); q_it != seqs.cend(); ++q_it) {
        cout << "\r" << 100*(q_it-seqs.cbegin())/count << "%";
        bool match = false;
        int rejects = 0;
        //bool false_negative = false;

        kmer_bits q_bitset(0);
        get_kmer_bitset(*q_it, q_bitset);

        for (auto c_it = centroids.cbegin(); 
                (c_it != centroids.cend()) && (rejects <= max_rejects); ++c_it) {
            
            size_t target_bits = (c_it->first).count();
            kmer_bits b = (c_it->first) & q_bitset;
            size_t set_bits = b.count();
            //cout << set_bits << " : " << (int) (target_bits/1.1) << endl;
            if (set_bits >= target_bits*(dist.threshold()-0.05)) {
                double d = dist.distance(*q_it, c_it->second);
                if (d >= dist.threshold()) {
                    match = true; // found cluster
                    break;
                }
                ++rejects;
            }

            //if (dist.distance(*q_it, c_it->second) >= dist.threshold())
            //    false_negative = true;

        }

        if (!match) {

            //*if (false_negative)
            //    ++neg_count;
            
            // add new centroid and write to stream in FASTA format
            centroids.push_back({q_bitset, *q_it});

            // write centroid entry to to clusters file
            fs_clusters << "C " << setw(6) << centroid_count++ << " "
                        << (*q_it).desc() << "\n";
            // write FASTA format to centroids file
            fs_centroids << ">" << (*q_it).desc() << '\n'
                         << (*q_it).to_string()   << '\n';

        }
    }

    //cout << endl << "False negatives: " << neg_count << endl;

    cout << endl;
    return centroid_count;
}

void Cluster::get_kmer_bitset(const Seq& s, bitset<KMER_BITSET>& b) {
    
    static const uint32_t k2 = 2 * KMER_LEN;
    static const uint32_t mask = pow(2, k2) - 1;
    
    for (size_t i = 0; i <= s.length() - KMER_LEN; ++i) {
        uint32_t kmer = 0;

        //memcpy(&kmer_l, (longer + (i/4)), 4);
        kmer = Distance::stream2int(s.data() + (i/4));
        kmer >>= (32 - 2*(i % 4) - k2);
        kmer &= mask;
        b.set(kmer);
    }
}