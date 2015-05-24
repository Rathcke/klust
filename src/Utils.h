#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>

#include "Seq.h"
#include "Distance.h"

namespace Utils {

/**
 * Return random int in interval [a,b].
 */
int get_rand(int a, int b);

/**
 * Given a DNA base ("ACGT"), return a random base different from the argument.
 */
char get_rand_base_not(char c);

void permute(std::vector<Seq>& seqs, int count, double ratio,
        std::ofstream& fs_cts);

void permute_chunks(std::vector<Seq>& seqs, int count, double ratio,
						std::ofstream& fs_cts, int chunk_size);

/**
 * Write a space separated distance matrix for the sequences in the given
 * vector to the given output stream, using the given Distance object.
 */
void print_matrix(std::vector<Seq>& seqs, std::ostream& fs_mat, Distance& dist);

/**
 * Given a vector of Seq, return a vector of Seq which contains sequences that
 * are mutually dissimilar, i.e. below the threshold similarity in the given
 * Distance object.
 */
std::vector<Seq> dissimilar_seqs(std::vector<Seq>& seqs, Distance& dist);

} // namespace Utils

#endif
