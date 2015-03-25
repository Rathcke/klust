#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <set>
#include <queue>

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
int Cluster::clust(fstream& fs_in, fstream& fs_centroids, fstream& fs_clusters,
        Distance& dist, int count, int max_rejects) {

    struct Centroid {
        struct Seq seq;
        set<string> kmers;  // set of occurring kmers in seq
    };

    vector<struct Centroid> centroids; // alternatively vector<pair<Seq, set<string>>>

    struct Seq s;
    int i = 0;

    while (IO::read_sequence(fs_in, s) && i < count) {
        bool match = false;
        ++i;

        set<string> query_kmers = dist.kmers(s); // kmers in query sequence

        priority_queue<pair<int, int>> q;  // non-increasingly sorted queue

        // calculate intersection and store cardinality and index in queue
        for (vector<struct Centroid>::const_iterator it = centroids.begin();
                it != centroids.end(); ++it) {
            vector<string> intersect;
            set_intersection(query_kmers.begin(), query_kmers.end(),
                             (it->kmers).begin(), (it->kmers).end(),
                             back_inserter(intersect));
            //cout << "Size of intersect: " << intersect.size() << endl;
            q.push({intersect.size(), it - centroids.begin()});
        }

        int rejects = 0;
        // loop through centroids
        while(!q.empty() && !match && rejects < max_rejects) {
            ++rejects;
            //cout << "q not empty" << endl;
            int ctr_index = (q.top()).second;
            Seq ctr = centroids[ctr_index].seq;
            q.pop();

            if (dist.compare(s.data, ctr.data)) {
                // write s belongs to centroids[i] to fs_clusters
                fs_clusters << ctr_index << ": " << s.data << endl;
                match = true; // found cluster
                break;
            }
        }

        if (!match) {
            // add new centroid and write to stream in FASTA format
            Centroid new_centroid = {s, dist.kmers(s)};
            centroids.push_back(new_centroid);

            fs_centroids << '>' << s.desc << endl
                         << s.data << endl;
        }
    }

    return centroids.size();
}
