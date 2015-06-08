#ifndef IO_H
#define IO_H

#include <string>
#include <fstream>
#include <list>
#include <vector>

#include "Cluster.h"
#include "Seq.h"

namespace IO {

/**
 * Read the specified number of sequences from the given filestream (in FASTA
 * format), construct a Seq object for each sequence and store in the given
 * vector. Return the number of read sequences.
 */
int read_seqs(std::ifstream& fs, std::vector<Seq>& seqs, int count);

int read_seqs2(std::istream &fs, std::vector<Seq>& seqs, int count);

/**
 * Read [DR]NA sequence from given filestream in given string.
 * Return true on success and false if there's no more to read.
 */
bool read_sequence(std::ifstream& fs, std::string& s);

/**
 * Print statistics about clustering result and memory usage.
 */
void print_stats(const std::vector<Seq>& seqs, const std::list<Centroid>& cts);

/**
 * Given a list of Centroid structs, write centroids in FASTA format
 * to the given file stream.
 */
void write_centroids(const std::list<Centroid>& cts,
        std::ofstream& fs_centroids);

/**
 * Given a list of Centroid structs, write clustering results to the given file
 * stream in a format inspired by UCLUST's .uc format.
 */
void write_clusters(const std::list<Centroid>& cts,
        std::ofstream& fs_clusters);

/**
 * Output HTML, to given file stream, containing Springy JavaScript code for
 * visualizing clustering result.
 */
void springy(const std::list<Centroid>& cts, std::ofstream& fs);

} // namespace IO

#endif
