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

Simply running klust with no arguments will print instructions for usage:

```
$ ./klust
Usage: ./klust <FASTA input file>

Options:
-o, --centroids file   Output FASTA file for centroids
-u, --clusters file    Output file for clustering results
-t, --id t             Set similarity threshold/identity to t (in [0,1])
-k k                   Set the k in k-mer, used in similarity metric
-c, --count n          Read and process n sequences
-d, --sort_decr        Sort sequences by decresing length
-i, --sort_incr        Sort sequences by increasing length
-l, --depth d          Set the depth of the tree of divides, i.e.
                       cluster 2^d subsets of sequences and combine
-m, --max_rejects x    Max number of rejects when searching for centroid
-s, --step_size s      Step s characters between k-mers when comparing

--springy file         Generate springy JavaScript code
```


### Example

The following example will cluster the sequences in the FASTA format file
`../data/SILVA_10k.fasta` which contains 10,000 RNA sequences (the first 10k
from the `SILVA_119_SSURef_tax_silva.fasta` dataset). The program is run with
default clustering parameters and the input sequences are sorted by increasing
sequence length. The centroids will be written to the file `centroids.fasta`
and the clustering results to the file `clusters`:

```sh
./klust ../data/SILVA_10k.fasta -o centroids.fasta -u clusters --sort_incr
```
