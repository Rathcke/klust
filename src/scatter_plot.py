#!/usr/bin/env python3.3

import numpy as np
from matplotlib import pyplot as plt

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

# read data
kmer = np.matrix(np.loadtxt('data/k6_distmat'))
lev = np.matrix(np.loadtxt('data/lev_distmat'))

# transform into array of upper tringle of matrix
kmer_indices = np.triu_indices_from(kmer)
kmer_a = np.asarray( kmer[kmer_indices] )[-1]
kmer_a = [x for x in kmer_a if x != 0]

lev_indices = np.triu_indices_from(lev)
lev_a = np.asarray( lev[lev_indices] )[-1]
lev_a = [x for x in lev_a if x != 0]

# scatter plot
plt.scatter(kmer_a, lev_a, s=20, c='g')

plt.show()
