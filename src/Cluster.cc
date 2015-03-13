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
int Cluster::clust(fstream& fs_in, fstream& fs_centroids,
        int threshold, int k, int count) {
    vector<struct Seq> centroids;
    struct Seq s;

    int i = 0;
    while (IO::readSequence(fs_in, s) && i < count) {
        ++i;
        if (matchCentroid(s, centroids, threshold, k)) {
            continue; // continue if s is close enough to some seq in centroids
        }
        centroids.push_back(s);
        fs_centroids << '>' << s.desc << endl
                     << s.data << endl;
    }
    return centroids.size();
}

int Cluster::clust(fstream& fs_in, fstream& fs_centroids, fstream& fs_clusters,
        int threshold, int k, int count) {
    vector<struct Seq> centroids; // TODO: gets big, maybe use more structure to speed up
    struct Seq s;

    int i = 0;
    while (IO::readSequence(fs_in, s) && i < count) {
        bool match = false;
        ++i;

        for (vector<Seq>::size_type i = 0; i < centroids.size(); ++i) {
            if (Distance::d2window(s.data, centroids[i].data, k, threshold)) {
                // write s belongs to centroids[i] to fs_clusters
                fs_clusters << s.data << " " << centroids[i].data << endl;
                match = true; // found cluster
                break;
            }
        }

        if (!match) {
            // add new centroid and write to stream in FASTA format
            centroids.push_back(s);
            fs_centroids << '>' << s.desc << endl
                         << s.data << endl;
        }
    }

    return centroids.size();
}

/**
 * Return true if given struct seq sequence, s, is within a given threshold, t,
 * of a cluster centroid in the given vector of struct seq, cs.
 */
bool Cluster::matchCentroid(const struct Seq& s, const vector<struct Seq>& cs,
                            int t, int k) {
    for (vector<struct Seq>::size_type i = 0; i < cs.size(); ++i) {
        if (Distance::d2window(s.data, cs[i].data, k, t)) {
            return true;
        }
    }
    return false;
}
