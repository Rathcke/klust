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
        Distance dist;
        IO io;

        bool matchCentroid(const struct seq& s,
                const std::vector<struct seq>& cs, int t, int k);
};

#endif
