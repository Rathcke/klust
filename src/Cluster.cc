#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <list>

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
    vector<Seq> centroids; // TODO: gets big, maybe use more structure to speed up
    map<key, list<pair<id, int>>> key_map;    
    struct Seq s;

    int i = 0;
    while (IO::read_sequence(fs_in, s) && i < count) {
        bool match = false;
        ++i;
        vector<pair<key,int>> s_keys = dist.compute_key(s.data, max_rejects);
        for (vector<pair<key,int>>::const_iterator it = s_keys.begin(); 
                it != s_keys.end(); ++it) {
            if (match) {
                break;
            }
            if (key_map.find(it->first) == key_map.end()) {
                continue;
            }
            list<pair<id,int>> lst = key_map.find(it->first)->second;
            for (list<pair<id,int>>::const_iterator it = lst.begin();
                    it != lst.end(); ++it) {
                if (dist.compare(s.data, centroids[it->first].data)) {
                    // write s belongs to centroids[i] to fs_clusters
                    fs_clusters << s.data << " " << centroids[it->first].data << endl;
                    match = true; // found cluster
                    break;
                }
            }
        }
    
        if (!match) {
            // add new centroid and write to stream in FASTA format
            vector<pair<key, int>> vec = dist.compute_key(s.data, max_rejects);
            if (vec.empty()) {
              throw logic_error("Calling compute_key on " 
                            + s.data + " returns an empty vector" );  
            }
            for (vector<pair<key,int>>::const_iterator it = vec.begin(); 
                    it != vec.end(); ++it) {
                int index = centroids.size();
                if (key_map.find(it->first) == key_map.end()) {
                    list<pair<id, int>> initlst;
                    key_map.insert({it->first, initlst});
                }
                list<pair<id,int>> ls = key_map.find(it->first)->second;
                list<pair<id,int>> elem;
                elem.push_back({index, it->second});
                ls.merge(elem); // TODO: Merges in increasing order - should be decreasing.
            }
            centroids.push_back(s);
            fs_centroids << '>' << s.desc << endl
                         << s.data << endl;
        }
    }
    for (map<key,list<pair<id,int>>>::const_iterator it = key_map.begin();
            it != key_map.end(); ++it) {
        int a = it->first;
        cout << "key: " << a << "  list: ";
        for (list<pair<id,int>>::const_iterator it2 = (it->second).begin();
                it2 != (it->second).end(); ++it2) {
            cout << it2->first << " ";
        }
        cout << endl;
    }
    return centroids.size();
}