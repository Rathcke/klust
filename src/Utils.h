#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>

#include "Seq.h"

namespace Utils {

int get_rand(int a, int b);

char get_rand_base_not(char c);

void permute(std::vector<Seq>& seqs, int count, double ratio,
        std::ofstream& fs_cts);

} // namespace Utils

#endif
