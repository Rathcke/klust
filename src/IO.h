#ifndef IO_H
#define IO_H

#include <string>
#include <fstream>

#include "Cluster.h"
#include "Seq.h"

namespace IO {

int read_seqs(std::ifstream& fs, std::vector<Seq>& seqs, int count);

bool read_sequence(std::ifstream& fs, std::string& s);

void write_results(const std::list<Centroid>& cts,
        std::ofstream& fs_centroids, std::ofstream& fs_clusters);

} // namespace IO

#endif
