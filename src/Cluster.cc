#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>

#include "Cluster.h"
#include "Distance.h"
#include "IO.h"

using namespace std;

/**
 * Given file stream, read sequences of FASTA format, compute the centroids of
 * the clustering (given threshold and k for k-mers) of the given count of
 * sequences. Output the centroids in FASTA format to the output file stream,
 * output clusters to file stream and return the number of centroids.
 */
int Cluster::simple_clust(fstream& fs_in, fstream& fs_centroids, fstream& fs_clusters,
        Distance& dist, int count, int max_rejects) {
    map<int, struct Seq> centroids; // TODO: gets big, maybe use more structure to speed up
    struct Seq s;

    int i = 0;
    int centroid_count = 0;

    while (IO::read_sequence(fs_in, s) && i < count) {
        bool match = false;
        ++i;
        vector<int> s_keys = dist.compute_key(s.data, max_rejects);

        for (vector<int>::const_iterator it = s_keys.begin(); 
                it != s_keys.end(); ++it) {
            if (centroids.find(*it) == centroids.end()) {
                continue;
            }
            if (dist.compare(s.data, (centroids.find(*it)->second).data)) {
                // write s belongs to centroids[i] to fs_clusters
                fs_clusters << s.data << " " << (centroids.find(*it)->second).data << endl;
                match = true; // found cluster
                break;
            }
        }

        if (!match) {
            // add new centroid and write to stream in FASTA format
            vector<int> vec = dist.compute_key(s.data, 1);
            if (vec.empty()) {
              throw logic_error("Calling compute_key on " 
                            + s.data + " returns an empty vector" );  
            }
            centroids.insert({vec[0], s});
            ++centroid_count;
            fs_centroids << '>' << s.desc << endl
                         << s.data << endl;
        }
    }

    return centroid_count;
}
