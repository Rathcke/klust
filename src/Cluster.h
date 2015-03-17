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

        static int clust(std::fstream& fs_in, std::fstream& fs_centroids,
                std::fstream& fs_clusters, Distance& dist, int count, 
                int max_rejects);

    private:
};

#endif
