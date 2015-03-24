#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <set>

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
    map<string, set<string>> centroids;
    struct Seq s;

    int i = 0;
    while (IO::read_sequence(fs_in, s) && i < count) {

        bool match = false;
        ++i;
        
        set<string> query_kmers = dist.kmers(s);
        for (set<string>::const_iterator it = query_kmers.begin();
            it != query_kmers.end(); ++it) {
            cout << *it << endl;
        }

        set<string> kmer_intersect;

        for (map<string, set<string>>::const_iterator it = centroids.begin();
                it != centroids.end(); ++it) {

           set_intersection(query_kmers.begin(), query_kmers.end(), 
            (it->second).begin(), (it->second).end(), kmer_intersect.begin());
        }


        if (!match) {
            // add new centroid and write to stream in FASTA format
            centroids.insert(pair<string, set<string>>(s.data, dist.kmers(s)));

            fs_centroids << '>' << s.desc << endl
                         << s.data << endl;
        }
    }
    return centroids.size();
}




/*
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
        }*/