#!/usr/bin/env python3

# Multidimensional scaling of synthetic dataset based on 380 sequences from
# SILVA which are at most 0.6 similar. 9 copies of each after the centroid.

import sys

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
count = 500
similarities = np.loadtxt('SILVA_nosort_500_distmat.txt')
#similarities = np.loadtxt('SILVA_sort_incr_500_distmat.txt')
#similarities = np.loadtxt('SILVA_sort_decr_500_distmat.txt')

# Perform Multi Dimensional Scaling
mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=0,
                   dissimilarity="precomputed", n_jobs=1)
pos = mds.fit(similarities).embedding_

# Perform isomap
#pos = manifold.Isomap(20, n_components=2).fit_transform(similarities)

# Perform t-SNE
tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
pos = tsne.fit_transform(similarities)

# Different color for each cluster
colors= ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6",
         "#A30059", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43",
         "#8FB0FF", "#997D87", "#5A0007", "#809693", "#1B4400", "#4FC601",
         "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900",
         "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#7B4F4B",
         "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", "#372101",
         "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
         "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68",
         "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
         "#BEC459", "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD",
         "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
         "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757",
         "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837",
         "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#201625",
         "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534",
         "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
         "#83AB58", "#001C1E", "#004B28", "#C8D0F6", "#A3A489", "#806C66",
         "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59",
         "#8ADBB4", "#1E0200", "#5B4E51", "#C895C5", "#320033", "#FF6832",
         "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58", "#7A7BFF",
         "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393",
         "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325",
         "#02525F", "#0AA3F7", "#E98176", "#DBD5DD", "#5EBCD1", "#3D4F44",
         "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5", "#E773CE",
         "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B",
         "#E704C4", "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7",
         "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058", "#A45B02",
         "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B",
         "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31",
         "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9",
         "#C6DC99", "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2",
         "#CCAA35", "#374527", "#8BB400", "#797868", "#C6005A", "#3B000A",
         "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
         "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6",
         "#71BB8C", "#00B433", "#789EC9", "#6D80BA", "#953F00", "#5EFF03",
         "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F", "#003109", "#0060CD",
         "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
         "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8",
         "#EA8B66", "#958A9F", "#BDC9D2", "#9FA064", "#BE4700", "#658188",
         "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71",
         "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
         "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46",
         "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C",
         "#92896B"]

cols = ['k'] * count

# first 100 of SILVA unsorted
#clusters = [(20, [64]), (0, []), (16, []), (66, []), (65, []),
#            (93, [52,54,55,59,51,56,9,58,10,11,60,68]), (40, []), (42, []),
#            (36, [34,5,35]), (85, []),
#            (88, [29,73,78,82,84,86,45,33,91,92,95,49,57,98,41,61,53,15,70,74,75]),
#            (31, []), (12, []), (47, [8]), (48, [50]), (18, [37,94]), (23, []),
#            (2, [22]), (28, [26]), (7, []), (67, [13]), (63, []), (76, []),
#            (72, []), (27, [32,3]), (14, []), (71, []), (62, []), (24, [19]),
#            (1, []), (39, []), (44, []), (46, []), (69, []), (4, []), (43, []),
#            (30, []), (6, [38]), (81, [77]), (89, []), (96, [25,79]),
#            (83, [17,80]), (99, [90,21,87]), (97, [])]

# first 500 of SILVA unsorted
clusters = [(499, []), (228, [229,314,329,370,422,423,424,425,472,487,491,493,496,498]), (497, []), (495, []), (420, [494]), (492, []), (411, [412,476,484,486,488,490]), (93, [157,160,166,258,273,286,354,384,429,448,459,489]), (421, [485]), (483, []), (413, [414,452,467,469,475,477,478,480,482]), (479, [481]), (202, [264,295,415,416,417,418,471,473,474]), (409, [470]), (408, [461,468]), (407, [466]), (206, [240,257,270,465]), (76, [249,375,385,400,403,404,405,406,450,451,453,454,460,462,464]), (463, []), (402, [458]), (457, []), (401, [456]), (455, []), (449, []), (447, []), (222, [223,224,225,226,293,302,305,315,318,324,326,361,364,446]), (445, []), (444, []), (443, []), (442, []), (440, [441]), (439, []), (438, []), (247, [437]), (436, []), (434, [435]), (333, [336,338,433]), (432, []), (431, []), (430, []), (299, [355,428]), (427, []), (426, []), (419, []), (410, []), (248, [372,399]), (374, [398]), (397, []), (363, [366,371,379,381,393,396]), (395, []), (394, []), (242, [392]), (391, []), (5, [34,35,36,390]), (237, [342,353,388,389]), (387, []), (231, [232,234,235,236,320,323,327,347,350,376,380,382,386]), (383, []), (233, [341,344,378]), (360, [377]), (97, [163,211,212,230,280,282,317,325,332,335,373]), (357, [369]), (368, []), (227, [308,311,367]), (348, [365]), (238, [239,345,359,362]), (221, [290,309,312,321,349,358]), (339, [356]), (296, [352]), (351, []), (346, []), (343, []), (340, []), (207, [208,209,213,214,215,216,259,267,269,271,274,275,276,277,279,285,288,291,294,297,304,307,316,319,322,331,334,337]), (26, [28,330]), (328, []), (261, [263,310,313]), (4, [306]), (217, [303]), (255, [301]), (300, []), (298, []), (252, [292]), (289, []), (287, []), (284, []), (39, [201,262,283]), (281, []), (278, []), (272, []), (204, [268]), (203, [266]), (265, []), (260, []), (137, [141,145,149,153,156,189,190,191,192,193,194,196,197,199,250,256]), (254, []), (253, []), (251, []), (241, [244,246]), (245, []), (243, []), (220, []), (219, []), (218, []), (210, []), (205, []), (200, []), (198, []), (195, []), (105, [109,113,117,121,129,133,159,161,162,164,167,168,169,170,171,172,174,175,177,178,179,180,181,182,183,184,185,187,188]), (186, []), (173, [176]), (165, []), (100, [103,146,154,158]), (77, [81,142,143,147,155]), (138, [139,140,148,152]), (85, [151]), (150, []), (69, [135,144]), (15, [29,33,41,45,49,53,57,61,70,73,74,75,78,82,84,86,88,91,92,95,98,104,106,107,115,116,118,119,122,123,124,126,128,130,132,134,136]), (131, []), (65, [127]), (125, []), (18, [37,94,102,108,110,111,112,120]), (114, []), (101, []), (87, [90,99]), (79, [96]), (89, []), (17, [21,80,83]), (72, []), (71, []), (9, [10,11,51,52,54,55,56,58,59,60,68]), (13, [67]), (66, []), (20, [64]), (63, []), (62, []), (48, [50]), (8, [47]), (46, []), (44, []), (43, []), (42, []), (40, []), (6, [38]), (3, [27,32]), (31, []), (30, []), (25, []), (19, [24]), (23, []), (2, [22]), (16, []), (14, []), (12, []), (7, []), (1, []), (0, [])]

# first 500 of SILVA incr. sorted
#clusters = [(498, [499]), (497, []), (496, []), (495, []), (494, []), (492, [493]), (98, [266,297,436,491]), (17, [18,19,22,49,76,187,188,199,223,317,451,458,479,480,481,482,483,484,485,486,487,488,489,490]), (442, [457,459,460,461,462,463,465,467,468,469,470,471,472,478]), (73, [74,80,88,97,132,144,155,220,230,248,268,304,305,335,352,477]), (476, []), (473, [475]), (100, [150,269,353,474]), (466, []), (464, []), (454, [456]), (448, [452,455]), (453, []), (314, [450]), (81, [449]), (447, []), (296, [438,439,440,443,446]), (445, []), (444, []), (311, [315,316,333,441]), (435, [437]), (140, [154,156,170,253,301,303,313,434]), (433, []), (359, [360,361,363,364,365,366,367,368,369,370,372,373,374,375,376,377,379,380,382,383,384,385,387,389,390,391,395,403,411,412,414,416,417,418,419,420,421,422,427,429,430,432]), (431, []), (428, []), (240, [261,426]), (102, [151,172,357,362,371,378,381,386,388,392,393,394,396,397,398,399,400,401,402,404,405,406,407,408,409,410,413,415,423,424,425]), (358, []), (356, []), (309, [355]), (354, []), (351, []), (347, [350]), (349, []), (348, []), (346, []), (16, [29,53,67,75,82,173,189,191,327,331,338,342,344,345]), (343, []), (340, [341]), (337, [339]), (193, [332,336]), (195, [321,322,324,325,326,328,329,334]), (152, [182,330]), (65, [87,95,96,108,111,323]), (320, []), (319, []), (35, [318]), (312, []), (310, []), (48, [57,58,59,60,62,66,272,302,308]), (307, []), (306, []), (298, [300]), (299, []), (169, [295]), (106, [258,270,290,294]), (107, [109,134,293]), (292, []), (263, [291]), (245, [274,275,276,277,278,280,281,282,283,284,285,286,287,288,289]), (192, [279]), (273, []), (271, []), (267, []), (265, []), (256, [264]), (136, [137,138,147,149,159,160,162,165,166,167,168,241,247,262]), (260, []), (259, []), (257, []), (246, [255]), (254, []), (231, [252]), (251, []), (250, []), (244, [249]), (204, [205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,224,225,226,227,228,229,235,236,239,242,243]), (0, [21,158,161,164,179,186,196,201,202,203,238]), (200, [237]), (234, []), (233, []), (1, [8,184,197,232]), (198, []), (110, [114,194]), (190, []), (185, []), (178, [180,181,183]), (177, []), (171, [175,176]), (174, []), (163, []), (157, []), (153, []), (139, [142,146,148]), (145, []), (143, []), (141, []), (135, []), (90, [91,92,127,133]), (130, [131]), (117, [129]), (85, [128]), (126, []), (113, [122,123,124,125]), (121, []), (112, [120]), (115, [119]), (118, []), (116, []), (86, [105]), (104, []), (77, [103]), (101, []), (99, []), (94, []), (93, []), (37, [40,89]), (84, []), (83, []), (23, [25,79]), (78, []), (72, []), (71, []), (70, []), (69, []), (68, []), (64, []), (63, []), (24, [61]), (54, [56]), (55, []), (52, []), (51, []), (50, []), (47, []), (46, []), (44, [45]), (43, []), (42, []), (38, [41]), (26, [27,30,31,34,39]), (36, []), (33, []), (32, []), (28, []), (20, []), (15, []), (14, []), (13, []), (5, [9,12]), (2, [3,4,6,7,10,11])]

# first 500 of SILVA decr. sorted
#clusters = [(478, [499]), (493, [494,496,498]), (497, []), (487, [495]), (269, [302,315,492]), (488, [490,491]), (489, []), (486, []), (485, []), (484, []), (483, []), (450, [477,480,481,482]), (479, []), (474, [476]), (475, []), (458, [465,468,469,472,473]), (471, []), (470, []), (467, []), (466, []), (464, []), (463, []), (410, [461,462]), (459, [460]), (457, []), (329, [347,456]), (455, []), (453, [454]), (452, []), (437, [438,439,441,444,451]), (449, []), (169, [448]), (447, []), (416, [432,446]), (443, [445]), (442, []), (440, []), (436, []), (435, []), (434, []), (191, [197,218,433]), (431, []), (430, []), (429, []), (428, []), (427, []), (371, [402,411,418,425,426]), (417, [424]), (284, [300,309,312,423]), (422, []), (421, []), (420, []), (419, []), (415, []), (374, [414]), (393, [413]), (172, [386,389,403,405,412]), (409, []), (373, [375,404,406,408]), (407, []), (401, []), (25, [146,230,349,400]), (399, []), (81, [84,87,92,95,96,97,98,99,100,101,102,103,104,105,107,108,109,110,111,112,113,114,115,135,136,140,141,142,154,155,158,162,168,330,398]), (224, [305,397]), (396, []), (380, [384,395]), (205, [210,229,240,394]), (307, [387,392]), (208, [370,390,391]), (379, [388]), (385, []), (383, []), (376, [382]), (381, []), (378, []), (377, []), (324, [365,367,369,372]), (366, [368]), (364, []), (336, [337,338,350,353,361,362,363]), (351, [352,357,360]), (359, []), (358, []), (356, []), (21, [148,165,194,196,231,250,266,285,346,355]), (354, []), (320, [348]), (345, []), (243, [327,342,344]), (264, [296,298,299,303,313,321,333,339,343]), (341, []), (340, []), (239, [251,265,328,332,334,335]), (204, [331]), (322, [325,326]), (323, []), (316, [317,318,319]), (314, []), (311, []), (173, [310]), (308, []), (163, [166,306]), (164, [170,171,174,175,176,177,304]), (301, []), (260, [297]), (257, [259,261,262,263,270,271,272,273,274,275,276,277,278,279,280,281,282,283,286,287,288,289,290,291,292,293,294,295]), (246, [268]), (267, []), (258, []), (66, [236,249,253,256]), (255, []), (245, [254]), (252, []), (248, []), (235, [247]), (244, []), (242, []), (241, []), (206, [238]), (237, []), (234, []), (8, [63,202,233]), (232, []), (228, []), (209, [211,212,213,214,215,216,217,219,220,221,223,225,226,227]), (222, []), (207, []), (203, []), (199, [201]), (200, []), (43, [47,52,53,56,59,60,61,65,186,195,198]), (193, []), (143, [192]), (190, []), (189, []), (188, []), (58, [167,181,185,187]), (184, []), (183, []), (179, [180,182]), (178, []), (160, [161]), (157, [159]), (156, []), (153, []), (152, []), (149, [151]), (150, []), (147, []), (145, []), (144, []), (75, [76,77,78,79,80,82,83,85,86,88,89,90,91,93,94,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,137,138,139]), (106, []), (67, [69,71,72,74]), (73, []), (70, []), (68, []), (62, [64]), (23, [26,28,31,32,33,34,35,36,37,38,39,40,42,57]), (55, []), (54, []), (51, []), (50, []), (9, [10,11,12,13,14,15,16,17,18,19,20,41,49]), (48, []), (44, [46]), (45, []), (30, []), (29, []), (24, [27]), (22, []), (6, [7]), (5, []), (4, []), (3, []), (2, []), (0, [1])]

#markers = [u'.'] * len(pos[:,0])
markers = [u'+'] * len(pos[:,0])
zorders = [1] * len(pos[:,0])
sizes = [40] * len(pos[:,0])
alphas = [0.9] * len(pos[:,0])

#for c in centroids:
#    #colors[c] = 'b'
#    markers[c] = u'.'
#    zorders[c] = 5
#    sizes[c] = 50
#    alphas[c] = 1.0

i = 0
for t, s in clusters:
    cols[t] = colors[i]
    markers[t] = u'.'
    zorders[t] = 5
    sizes[t] = 50
    alphas[t] = 1.0

    for e in s:
        cols[e] = colors[i]
    i += 1

#plt.figure(num=None, figsize=(16, 12), dpi=100, facecolor='w', edgecolor='k')

# Scatter plot
#plt.scatter(pos[:,0], pos[:,1], s=20, c=list(colors))
for x, y, color, m, sz, z, a in zip(pos[:,0], pos[:,1], cols, markers,
        sizes, zorders, alphas):
    plt.scatter(x, y, s=sz, c=color, marker=m, zorder=z, alpha=a)

# Label centroids
#c = 0
#for x, y in zip(pos[:,0], pos[:,1]):
#    if c in centroids:
#        plt.annotate(
#            str(c),
#            xy= (x,y), xytext=(-20,20),
#            textcoords='offset points', ha='right', va='bottom',
#            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0') )
#    c += 1

#plt.axis('off')

plt.xticks([])
plt.yticks([])

# no sort
plt.xlim(-28.0, 36.0)
plt.ylim(-35.0, 38.0)

# sort_decr
#plt.xlim(-30.0, 36.0)
#plt.ylim(-35.0, 40.0)

# sort_incr
#plt.xlim(-30.0, 38.0)
#plt.ylim(-35.0, 40.0)

plt.tight_layout()

plt.savefig("SILVA_t-SNE_" + str(count) + ".svg",
        bbox_inches="tight")

#plt.show()
