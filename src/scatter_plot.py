#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

# read data
kmer = np.matrix(np.loadtxt('data/k8_distmat'))
lev = np.matrix(np.loadtxt('data/lev_distmat'))

# transform into array of upper tringle of matrix
kmer_indices = np.triu_indices_from(kmer)
kmer_a = np.asarray( kmer[kmer_indices] )[-1]
kmer_a = [x for x in kmer_a if x != 0]

lev_indices = np.triu_indices_from(lev)
lev_a = np.asarray( lev[lev_indices] )[-1]
lev_a = [x for x in lev_a if x != 0]

# scatter plot
#plt.scatter(kmer_a, lev_a, s=20, c='g')
plt.scatter(kmer_a, lev_a, s=1, c='g', alpha=0.1)

#plt.suptitle('test title', fontsize=20)
plt.xlabel('k-mer similarity (k = 8)', fontsize=18)
plt.ylabel('Levenshtein similarity', fontsize=18)

plt.savefig("k8.png")
#plt.show()
