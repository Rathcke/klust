#ifndef CLUSTER_H
#define CLUSTER_H

#include <string>
#include <vector>
#include <bitset>
#include <list>

#include <functional>

#include "Distance.h"

#define KMER_BITSET 4096
#define KMER_LEN 6

struct Centroid {
    Seq& seq;                       // centroid sequence
    std::bitset<KMER_BITSET> bits;  // bitset of the k-mers occurring in seq
    size_t count;                   // # of distinct k-mers in seq

    // sequences in the cluster represented by the centroid
    std::vector<std::reference_wrapper<Seq>> cls_seqs;

    const unsigned num; // numbering of centroid in the order of discovery

    Seq *link;

    Centroid(Seq& s, std::bitset<KMER_BITSET> bits, unsigned int num)
            : seq {s}, bits {bits}, num {num} {
        count = bits.count(); // # of distinct kmers in seq
        link = nullptr;
    }
};

class Cluster
{
    public:
        Cluster(Distance&& d, int max_rejects)
            : dist {d}, max_rejects {max_rejects} {}

        int simple_clust(const std::vector<Seq>& seqs, std::ofstream& fs_centroids,
                std::ofstream& fs_clusters, Distance& dist, int count,
                int max_rejects);

        int thorough_clust(const std::vector<Seq>& seqs,
                std::ofstream& fs_centroids, std::ofstream& fs_clusters,
                Distance& dist, int count);

        void kmer_select_clust(std::vector<Seq>::iterator begin,
                std::vector<Seq>::iterator end, std::list<Centroid>& cts);

        int clust(std::vector<Seq>& seqs, std::list<Centroid>& cts, int depth = 0);

    private:
        Distance& dist;
        int max_rejects;

        void merge(std::list<Centroid>& res, const std::list<Centroid>& c1);
};

#endif
