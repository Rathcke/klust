#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <list>
#include <sys/time.h>

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
    typedef int id;
    typedef int key;
    vector<Seq> centroids;
    map<key, list<pair<id, int>>> key_map;    
    struct Seq s;

    int i = 0;
    while (IO::read_sequence(fs_in, s) && i < count) {
        bool match = false;
        ++i;
        vector<pair<key,int>> s_keys = dist.compute_key(s.data, max_rejects);
        if (s_keys.empty()) {
          throw logic_error("Calling compute_key on " 
                        + s.data + " returns an empty vector");  
        }
        int rejects = 0;

        // loop through most common k-mers (i.e. lexicographically ordered
        // keys for k-mers) in s:
        for (vector<pair<key,int>>::const_iterator it = s_keys.begin(); 
                it != s_keys.end(); ++it) {
            if (match || rejects >= max_rejects) {    // we already found a matching centroid
                break;
            }
            if (key_map.find(it->first) == key_map.end()) {
                continue;   // the k-mer doesn't exists in the centroid map
            }
            
            // loop through list of sequences with this k-mer as a common key;
            // stop when finding a matching centroid or after max_rejects tries
            list<pair<id,int>> lst = key_map.find(it->first)->second;
            for (list<pair<id,int>>::const_iterator it = lst.begin();
                    it != lst.end() && rejects < max_rejects; ++it, ++rejects) {
                if (dist.compare(s.data, centroids[it->first].data)) {
                    // write s belongs to centroids[i] to fs_clusters
                    fs_clusters << it->first << ": " << s.data << endl;
                    match = true; // found cluster
                    break;
                }

            }
        }              

        if (!match) {
            // add new centroid and write to stream in FASTA format
            for (vector<pair<key,int>>::const_iterator it = s_keys.begin(); 
                    it != s_keys.end(); ++it) {
                int index = centroids.size();
                if (key_map.find(it->first) == key_map.end()) {
                    // initial list with 1 times (index, key_count) tuple
                    list<pair<id, int>> initlst(1, {index,it->second});
                    key_map.insert({it->first, initlst});
                } else {
                    list<pair<id,int>> elem;
                    elem.push_back({index, it->second});
                    (key_map.find(it->first)->second).merge(elem,
                        [](const pair<id, int>& lhs, const pair<id, int>& rhs) {
                            return lhs.second > rhs.second;
                        });
                }
            }
            centroids.push_back(s);
            fs_centroids << '>' << s.desc << endl
                         << s.data << endl;
        }
    }

/*    for (map<key,list<pair<id,int>>>::const_iterator it = key_map.begin();
            it != key_map.end(); ++it) {
        int a = it->first;
        cout << "key: " << a << "  list: ";
        for (list<pair<id,int>>::const_iterator it2 = (it->second).begin();
                it2 != (it->second).end(); ++it2) {
            cout << "(" << it2->first << ", " << it2->second << ") ";
        }
        cout << endl;
    }*/
    return centroids.size();
}
