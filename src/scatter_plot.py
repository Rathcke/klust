#!/usr/bin/env python3.3

import numpy as np
from matplotlib import pyplot as plt

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

# read data
#k4  = np.loadtxt('data/k4_distmat')
#lev = np.loadtxt('data/lev_distmat')
k4  = np.matrix(np.loadtxt('data/k4_distmat'))
lev = np.matrix(np.loadtxt('data/lev_distmat'))

k4_indices = np.triu_indices_from(k4)
k4_a = np.asarray( k4[k4_indices] )[-1]

lev_indices = np.triu_indices_from(lev)
lev_a = np.asarray( lev[k4_indices] )[-1]

#tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
#pos = tsne.fit_transform(similarities)

#Scatter plot
#plt.scatter(pos[:, 0], pos[:, 1], s=20, c='g')
plt.scatter(k4_a, lev_a, s=20, c='g')

plt.show()
