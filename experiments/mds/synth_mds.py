#!/usr/bin/env python3

# Multidimensional scaling of synthetic dataset based on 380 sequences from
# SILVA which are at most 0.6 similar. 9 copies of each after the centroid.

import numpy as np

from itertools import cycle

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.cm as cmx
import matplotlib.colors as colors

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

# Read data
count = 1600
similarities_all = np.loadtxt('distmat_synth_SILVA_3800_dist.txt')
similarities = similarities_all[0:count,0:count]

# Perform Multi Dimensional Scaling
mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=0,
                   dissimilarity="precomputed", n_jobs=1)
pos = mds.fit(similarities).embedding_

# Perform isomap
#pos = manifold.Isomap(20, n_components=2).fit_transform(similarities)

# Perform t-SNE
tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
pos = tsne.fit_transform(similarities)

# Centroids every 10'th sequence
#centroids = range(0, 3800, 10)
centroids = range(0, count, 10)

# different color for each cluster
colors = 10 * 'b' + 10 * 'g' + 10 * 'r' + 10 * 'c' + 10 * 'm' + 10 * 'y' + 10 * 'k'
colors = cycle(colors)

markers = [u'.'] * len(pos[:,0])
zorders = [1] * len(pos[:,0])
sizes = [40] * len(pos[:,0])
alphas = [0.5] * len(pos[:,0])

for c in centroids:
    #colors[c] = 'r'
    markers[c] = u'+'
    zorders[c] = 5
    sizes[c] = 50
    alphas[c] = 1.0

# Scatter plot
#plt.scatter(pos[:,0], pos[:,1], s=20, c=list(colors))

for x, y, m, sz, z in zip(pos[:,0], pos[:,1], markers, sizes, zorders):
    plt.scatter(x, y, s=sz, c=next(colors), marker=m, zorder=z)

#plt.savefig("MDS_t-SNE_SILVA_500.png")
plt.show()
