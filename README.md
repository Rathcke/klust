# klust

Efficient DNA/RNA sequence clustering using *k*-mers as an approximation for
sequence similarity.


## Compiling and running klust

The klust program has been developed and tested on GNU/Linux systems and
compiled with the `g++` compiler (version 5.1.0), which is part of GCC. The
program has also been tested to work with the `clang++` compiler (version
3.6.1). The program is expected to compile and run on any POSIX compliant
system with a recent `g++` or `clang++` compiler supporting C++11.

The following commands will download and compile the `klust` program, assuming
`git` is installed:

```sh
git clone https://github.com/Rathcke/klust.git
cd klust/src
make
```


## Example

The following example will cluster the sequences in the FASTA format file
`../data/SILVA_10k.fasta` which contains 10,000 RNA sequences (the first 10k
from the `SILVA_119_SSURef_tax_silva.fasta` dataset). The program is run with
default clustering parameters and the input sequences are sorted by increasing
sequence length. The centroids will be written to the file `centroids.fasta`
and the clustering results to the file `clusters`:

```sh
./klust ../data/SILVA_10k.fasta centroids.fasta clusters --sort_incr
```
