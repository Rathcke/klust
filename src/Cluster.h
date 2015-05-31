#ifndef CLUSTER_H
#define CLUSTER_H

#include <bitset>
#include <functional>
#include <list>
#include <string>
#include <vector>

#include "Distance.h"

#define KMER_BITSET 4096
#define KMER_LEN 6

struct Centroid {
    const Seq& seq;                 // centroid sequence
    std::bitset<KMER_BITSET> bits;  // bitset of the k-mers occurring in seq
    size_t count;                   // # of distinct k-mers in seq

    // sequences in the cluster represented by the centroid
    std::vector<std::reference_wrapper<const Seq>> cls_seqs;

    const unsigned num; // numbering of centroid in the order of discovery

    Centroid *link;

    Centroid(const Seq& s, std::bitset<KMER_BITSET> bits, unsigned int num)
            : seq {s}, bits {bits}, num {num} {
        count = bits.count(); // # of distinct kmers in seq
        link = nullptr;
    }
};

class Cluster
{
    public:
        Cluster(Distance& d, int max_rejects)
            : dist {d}, max_rejects {max_rejects} {}

        /**
         * Cluster the given vector of sequences and output the centroids in
         * FASTA format to the first file stream and the clustering results to
         * the second file stream. Return the number of clusters. Clustering is
         * done using the greedy Simple-Clust algorithm, which compares each
         * sequence with centroids with the same most frequently occurring
         * k-mers as the sequence.
         */
        int simple_clust(const std::vector<Seq>& seqs,
                std::ofstream& fs_centroids, std::ofstream& fs_clusters);

        /**
         * Cluster the sequences in the given collection using a naive
         * clustering algorithm, which compares each sequence sequence with
         * every centroid until a match is found. If not match is found, the
         * sequence becomes a new centroid. The centroids are written in FASTA
         * format to the first file stream, the clustering results are written
         * to the second file stream and the number of centroids is returned.
         */
        int thorough_clust(const std::vector<Seq>& seqs,
                std::ofstream& fs_centroids, std::ofstream& fs_clusters);

        /**
         * For every sequence in the given collection, search through the
         * centroids for one where the number of distinct k-mers in both the
         * query sequence and centroid are at least dist.threshold() times the
         * number of distinct k-mers in the centroid sequence. If no match is
         * found after max_rejects tries, the sequence becomes a new centroid.
         * The list of centroids is returned in the given argument.
         */
        void kmer_select_clust(std::vector<Seq>::const_iterator begin,
                std::vector<Seq>::const_iterator end, std::list<Centroid>& cts);

        /**
         * Given a collection of sequences, cluster 2^depth subparts of the
         * sequences in separate threads using the kmer_select_clust function,
         * combine the resulting list of centroids and return in the Centroid
         * list given as argument.
         */
        int clust(const std::vector<Seq>& seqs, std::list<Centroid>& cts, int depth = 0);

    private:
        Distance& dist;
        int max_rejects;

        /**
         * Merge two vectors of centroids and store the result in the first vector.
         */
        void merge(std::list<Centroid>& res, const std::list<Centroid>& c1);

        /**
         * Given a sequence, return the uint32_t representations of up to the n
         * most frequent k-mers or as many as exists.
         */
        std::vector<uint32_t> most_frequent_kmers(const Seq& s, int n);
};

#endif
