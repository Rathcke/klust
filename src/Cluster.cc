#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Cluster.h"
#include "Distance.h"
#include "IO.h"

using namespace std;

/**
 * Given file stream, read sequences of FASTA format, compute the centroids of
 * the clustering (given threshold and k for k-mers) of the given count of
 * sequences. Output the centroids in FASTA format to the output file stream
 * and return the number of centroids.
 */
int Cluster::clust(fstream& in, fstream& out, int threshold, int k, int count) {
    vector<struct seq> centroids;
    struct seq s;

    int i = 0;
    while (io.readSequence(in, s) && i < count) {
        ++i;
        if (matchCentroid(s, centroids, threshold, k)) {
            continue; // continue if s is close enough to some seq in centroids
        }
        centroids.push_back(s);
        out << '>' << s.desc << endl
            << s.data << endl;
    }
    return centroids.size();
}

/**
 * Return true if given struct seq sequence, s, is within a given threshold, t,
 * of a cluster centroid in the given vector of struct seq, cs.
 */
bool Cluster::matchCentroid(const struct seq& s, const vector<struct seq>& cs,
                            int t, int k) {
    for (vector<struct seq>::size_type i = 0; i < cs.size(); ++i) {
        if (dist.d2window(s.data, cs[i].data, k) <= t) {
            return true;
        }
    }
    return false;
}
