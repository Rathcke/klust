#ifndef CLUSTER_H
#define CLUSTER_H

#include <string>
#include <vector>
#include <bitset>

#include <functional>

#include "Distance.h"

#define KMER_BITSET 4096
#define KMER_LEN 6

struct Centroid {
    Seq& seq;
    std::bitset<KMER_BITSET> bits;
    size_t count;
    std::vector<std::reference_wrapper<Seq>> cls_seqs;
    //std::vector<Seq> cls_seqs; TODO: maybe move?

    Centroid(Seq& s, std::bitset<KMER_BITSET> bits) : seq(s), bits(bits) {
        count = bits.count(); // # of distinct kmers in seq
    }
};

class Cluster
{
    public:
        Cluster(Distance& d, int max_rejects)
            : dist {d}, max_rejects {max_rejects} {}

        /*static int intersect_clust(std::fstream& fs_in, std::fstream& fs_centroids,
                std::fstream& fs_clusters, Distance& dist, int count,
                int max_rejects);*/

        static int simple_clust(const std::vector<Seq>& seqs, std::ofstream& fs_centroids,
                std::ofstream& fs_clusters, Distance& dist, int count,
                int max_rejects);

        static int thorough_clust(const std::vector<Seq>& seqs,
                std::ofstream& fs_centroids, std::ofstream& fs_clusters,
                Distance& dist, int count);

        void kmer_select_clust(std::vector<Seq>::iterator begin,
                std::vector<Seq>::iterator end, std::vector<Centroid>& cts);

        //int clust(std::vector<Seq>& seqs, int subclusterings);

        int clust(std::vector<Seq>::iterator begin, std::vector<Seq>::iterator end,
                std::vector<Centroid>& cts, int depth);

        int kmer_clust(std::vector<Seq>& seqs, std::vector<Centroid>& cts);

    private:
        Distance& dist;
        int max_rejects;

        void merge(std::vector<Centroid>& c0, const std::vector<Centroid>& c1);
};

#endif
