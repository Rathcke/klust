#!/usr/bin/env python3

import sys

import numpy as np
from matplotlib import pyplot as plt

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

# read data
i, x, y = np.loadtxt(sys.argv[1])

# scatter plot
plt.scatter(i, x, s=30, c='r', alpha=0.6, marker=u'o')
plt.scatter(i, y, s=30, c='b', alpha=0.8, marker=u'+')

plt.savefig("cluster_counts.png")
#plt.show()
