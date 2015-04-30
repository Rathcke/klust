#ifndef IO_H
#define IO_H

#include <string>
#include <fstream>

#include "Seq.h"

namespace IO {

int read_seqs(std::ifstream& fs, std::vector<Seq>& seqs, int count);

bool read_sequence(std::ifstream& fs, std::string& s);

} // namespace IO

#endif
