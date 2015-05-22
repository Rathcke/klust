#ifndef IO_H
#define IO_H

#include <string>
#include <fstream>
#include <list>
#include <vector>

#include "Cluster.h"
#include "Seq.h"

namespace IO {

int read_seqs2(std::istream &fs, std::vector<Seq>& seqs, int count);

int read_seqs(std::ifstream& fs, std::vector<Seq>& seqs, int count);

bool read_sequence(std::ifstream& fs, std::string& s);

void print_stats(const std::vector<Seq>& seqs, const std::list<Centroid>& cts);

void write_centroids(const std::list<Centroid>& cts, std::ofstream& fs_centroids);

void write_clusters(const std::list<Centroid>& cts, std::ofstream& fs_clusters);

void springy(const std::list<Centroid>& cts, std::ofstream& fs);

} // namespace IO

#endif
