#ifndef IO_H
#define IO_H

#include <string>

class IO
{
    public:
        IO() {}

        bool readSequence(std::fstream& fs, std::string& s);

    private:

};

#endif
