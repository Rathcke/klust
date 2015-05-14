import numpy as np

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

# Read data
similarities = np.loadtxt('distmat.txt')

# Perform Multi Dimensional Scaling
mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=0,
                   dissimilarity="precomputed", n_jobs=1)
pos = mds.fit(similarities).embedding_

# Perform isomap
#pos = manifold.Isomap(20, n_components=2).fit_transform(similarities)

# Perform t-SNE
tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
pos = tsne.fit_transform(similarities)

# Scatter plot
plt.scatter(pos[:, 0], pos[:, 1], s=20, c='g')

# Centroids for t = 0.85
centroids = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        23, 25, 26, 29, 30, 31, 37, 39, 40, 42, 43, 44, 46, 48, 50, 62, 63, 65,
        66, 69, 71, 72, 76, 77, 79, 85, 87, 89, 93, 96, 97, 100, 101, 105, 110,
        114, 125, 131, 137, 141, 145, 150, 153, 165, 173, 176, 186, 195, 200,
        201, 202, 203, 204, 205, 206, 207, 210, 217, 218, 219, 220, 221, 222,
        227, 228, 231, 233, 235, 237, 238, 241, 242, 243, 245, 247, 248, 249,
        251, 252, 253, 254, 255, 260, 261, 264, 265, 272, 273, 278, 281, 283,
        284, 287, 289, 296, 298, 299, 300, 306, 328, 330, 333, 338, 339, 340,
        342, 343, 346, 347, 348, 351, 357, 360, 363, 368, 369, 374, 383, 387,
        388, 391, 394, 395, 397, 401, 402, 407, 408, 409, 410, 411, 413, 419,
        420, 421, 422, 426, 427, 430, 431, 432, 433, 434, 436, 438, 439, 440,
        442, 443, 444, 445, 447, 449, 452, 455, 457, 463, 475, 479, 483, 487,
        492, 493, 495, 497, 499]

# Centroids for t = 0.6
#centroids = [0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 14, 15, 17, 19, 20, 23, 30, 44, 46,
#        62, 65, 66, 71, 77, 97, 114, 131, 195, 202, 203, 220, 221, 227, 243,
#        253, 261, 287, 293, 296, 299, 351, 394, 407, 430, 432, 439, 445, 449]

# Centroids for thorough_clust and t = 0.6
#centroids = [0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 14, 15, 17, 19, 20, 23, 30, 44, 46,
#        62, 65, 66, 71, 77, 97, 114, 131, 195, 202, 203, 220, 221, 227, 243,
#        261, 287, 293, 296, 299, 351, 394, 407, 430, 432, 439, 445, 449, 463]

# Label centroids
c = 0
for x, y in zip(pos[:,0], pos[:,1]):
    if c in centroids:
        plt.annotate(
            str(c),
            xy= (x,y), xytext=(-20,20),
            textcoords='offset points', ha='right', va='bottom',
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0') )
    c += 1

plt.savefig("MDS_t-SNE_SILVA_500.png")
#plt.show()
