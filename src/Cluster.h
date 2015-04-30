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
    //std::vector<Seq> cls_seqs;    // TODO: maybe move?

    Centroid(Seq& s, std::bitset<KMER_BITSET> bits) : seq(s), bits(bits) {
        count = bits.count(); // # of distinct kmers in seq
    }

    /*Centroid& operator=(Centroid& other) {
        seq = std::ref(other.seq);
        bits = other.bits;
        count = other.count;
        cls_seqs = std::ref(other.cls_seqs);
        return *this;
    }*/
};

class Cluster
{
    public:
        Cluster() {}

        /*static int intersect_clust(std::fstream& fs_in, std::fstream& fs_centroids,
                std::fstream& fs_clusters, Distance& dist, int count,
                int max_rejects);*/

        static int simple_clust(const std::vector<Seq>& seqs, std::ofstream& fs_centroids,
                std::ofstream& fs_clusters, Distance& dist, int count,
                int max_rejects);

        static int thorough_clust(const std::vector<Seq>& seqs,
                std::ofstream& fs_centroids, std::ofstream& fs_clusters,
                Distance& dist, int count);

        static int kmers_select_clust(const std::vector<Seq>& seqs, std::ofstream& fs_centroids,
            std::ofstream& fs_clusters, Distance& dist, int max_rejects);


        static int clust(std::vector<Seq>& seqs, Distance& dist, int max_rejects);

    private:
        /*static int sub_clust(std::vector<Seq>::const_iterator begin,
                std::vector<Seq>::const_iterator end,
                std::vector<Centroid>& cts, Distance& dist, int max_rejects);*/

};

#endif
