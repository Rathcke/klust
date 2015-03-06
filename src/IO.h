#ifndef IO_H
#define IO_H

#include <string>
#include <fstream>

struct seq {
    std::string desc; // description following '>' in fasta format
    std::string data; // actual sequence data
};

class IO
{
    public:
        IO() {}

        bool readSequence(std::fstream& fs, std::string& s);
        bool readSequence(std::fstream& fs, struct seq& s);

    private:

};

#endif
