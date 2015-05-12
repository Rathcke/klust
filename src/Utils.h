#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>

#include "Seq.h"

namespace Utils {

void permute(std::vector<Seq>& seqs, int count, double ratio,
        std::ofstream& fs_cts);

} // namespace Utils

#endif
