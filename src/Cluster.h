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

        /*static int clust(std::fstream& fs_in, std::fstream& fs_centroids,
                int threshold, int k, int count);*/

        static int clust(std::fstream& fs_in, std::fstream& fs_centroids,
                std::fstream& fs_clusters, Distance& dist, int count);

    private:
        /*static bool matchCentroid(const struct Seq& s,
                const std::vector<struct Seq>& cs, int t, int k);*/
};

#endif
