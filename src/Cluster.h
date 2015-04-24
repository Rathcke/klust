#ifndef CLUSTER_H
#define CLUSTER_H

#include <string>
#include <vector>
#include <bitset>

#include "IO.h"
#include "Distance.h"

#define KMER_BITSET 4096
#define KMER_LEN 6

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

    private:

        static void get_kmer_bitset(const Seq& s, std::bitset<KMER_BITSET>& b);
};

#endif
