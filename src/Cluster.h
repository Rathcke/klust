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
    const Seq& seq;                 // centroid sequence

    // sequences in the cluster represented by the centroid
    std::vector<std::reference_wrapper<const Seq>> cls_seqs;

    const unsigned num; // numbering of centroid in the order of discovery

    const Seq *link;

    Centroid(const Seq& s, unsigned int num)
            : seq {s}, num {num} {
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

        void kmer_select_clust(std::vector<Seq>::const_iterator begin,
                std::vector<Seq>::const_iterator end, std::list<Centroid>& cts);

        int clust(const std::vector<Seq>& seqs, std::list<Centroid>& cts, int depth = 0);

    private:
        Distance& dist;
        int max_rejects;

        void merge(std::list<Centroid>& res, const std::list<Centroid>& c1);
};

#endif
