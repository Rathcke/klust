#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

# read data
kmer = np.matrix(np.loadtxt('../../data/k8_distmat'))
lev = np.matrix(np.loadtxt('../../data/lev_distmat'))

# transform into array of upper tringle of matrix
kmer_indices = np.triu_indices_from(kmer)
kmer_a = np.asarray( kmer[kmer_indices] )[-1]
#kmer_a = [x for x in kmer_a if x != 0]

lev_indices = np.triu_indices_from(lev)
lev_a = np.asarray( lev[lev_indices] )[-1]
#lev_a = [x for x in lev_a if x != 0]

# scatter plot
#plt.scatter(kmer_a, lev_a, s=20, c='g')
#plt.scatter(kmer_a, lev_a, s=1, c='g', alpha=0.1)
#plt.scatter(kmer_a, lev_a, s=2, c='g', alpha=0.4)

fig, ax = plt.subplots() #figsize=[5,4])

ax.scatter(kmer_a, lev_a, s=2, c='g', alpha=0.4)

ax.set_xlim(0.00, 1.05)
ax.set_ylim(0.02, 1.05)

#plt.suptitle('test title', fontsize=20)
plt.xlabel('k-mer similarity (k = 8)', fontsize=18)
plt.ylabel('Levenshtein similarity', fontsize=18)

# subplot
axins = zoomed_inset_axes(ax, 2.5, loc=4)
axins.scatter(kmer_a, lev_a, s=2, c='g', alpha=0.8)
axins.set_xlim(0.80,1.01)
axins.set_ylim(0.80,1.01)
axins.set_xticks([0.85, 0.90, 0.95, 1.0])
axins.set_yticks([0.85, 0.90, 0.95, 1.0])

mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5")

#axins.axis([0.85,1.01,0.85,1.01])
axins.xaxis.tick_top()

#m, b = np.polyfit(kmer_a, lev_a, 1)
#ax.plot(kmer_a, m*kmer_a + b, '-')

plt.draw()
plt.show()

#plt.savefig("distance_comparison_k8.png")
#plt.show()
