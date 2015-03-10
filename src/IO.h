#ifndef IO_H
#define IO_H

#include <string>
#include <fstream>

struct Seq {
    std::string desc; // description following '>' in fasta format
    std::string data; // actual sequence data
};

class IO
{
    public:
        IO() {}

        static bool readSequence(std::fstream& fs, std::string& s);
        static bool readSequence(std::fstream& fs, struct Seq& s);

    private:

};

#endif
