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

        static int clust(std::fstream& in, std::fstream& out,
                int threshold, int k, int count);

    private:
        static bool matchCentroid(const struct Seq& s,
                const std::vector<struct Seq>& cs, int t, int k);
};

#endif
