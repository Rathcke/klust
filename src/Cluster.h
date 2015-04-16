#ifndef CLUSTER_H
#define CLUSTER_H

#include <string>
#include <vector>

#include "IO.h"
#include "Distance.h"

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

    private:
};

#endif
