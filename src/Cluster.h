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

        int clust(std::fstream& in, std::fstream& out, int threshold, int k, int count);

    private:
        bool matchCentroid(std::string s, 
                const std::vector<std::string>& cs, int t, int k);
        IO io;
        Distance dist;
};

#endif
