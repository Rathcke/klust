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


int Cluster::clust(fstream& in, fstream& out, int threshold, int k, int count) {
    vector<string> centroids;
    string tmp;
    int i = 0;
    while (io.readSequence(in, tmp) && i < count) {
        i++;
        if (matchCentroid(tmp, centroids, threshold, k)) {
            continue;
        }
        centroids.push_back(tmp);
    }
    return centroids.size();
}


bool Cluster::matchCentroid(string s, const vector<string>& cs, int t, int k) {
    for (unsigned int i = 0; i < cs.size(); i++) {
        if (dist.d2window(s, cs[i], k) <= t) {
            return true;
        }
    }
    return false;
}