
# import general modules
import scprep
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import phate
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import csv
import math
import matplotlib.colors as colors

# import my functions
from ysc_Gresham2019.functions.getValuesForGeneGroups import getValuesForGeneGroups
from ysc_Gresham2019.functions.addORForGene import addORForGene
from ysc_Gresham2019.functions.getLogFcExpression import getLogFcExpression
from ysc_Gresham2019.functions.getDiffExpressedGenes import getDiffExpressedGenes
# import my scripts
from ysc_Gresham2019.scripts.clusteredCorrMatrix import clusteredCorrMatrix
from ysc_Gresham2019.scripts.getKmeansClusteringOnSelectedFeatures import getKmeansClusteringOnSelectedFeatures
from ysc_Gresham2019.scripts.externalDataAnalysis import externalDataAnalysis
from ysc_Gresham2019.scripts.similarityToExternalData import similarityToExternalData

# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

############################################################################################################

# 1. For MAGIC data: make correlation plot for all genes and get the list of gene clusters
YPD_magic_t2 = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2.csv', index_col=0)
[YPD_magic_t2__SlowRearranged, YPD_magic_t2_corrMtrxRearranged, YPD_magic_t2_myClusteredGroups] = \
    clusteredCorrMatrix(YPD_magic_t2, 'corr_AllGenes_MAGIC.png', 'MAGIC, all genes', threshGrowthRate=0.2)
# save rearranged dataset and clusters themselves
YPD_magic_t2__SlowRearranged.to_csv('./ysc_Gresham2019/Data/YPD_magic_t2__SlowRearranged.csv')
for i in range(len(YPD_magic_t2_myClusteredGroups)):
    with open('./ysc_Gresham2019/Data/MAGIC_clusteredGeneSets/geneSet_' + str(i) + '.csv', 'w', newline='\n') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(YPD_magic_t2_myClusteredGroups[i])
allGroupedGenes = YPD_magic_t2__SlowRearranged.columns
with open('./ysc_Gresham2019/Data/MAGIC_clusteredGeneSets/allGenes.csv', 'w', newline='\n') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(allGroupedGenes)

YPD_magic_t2__SlowRearranged = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2__SlowRearranged.csv', index_col=0)
stat_YPD_magic_t2__SlowRearranged = getValuesForGeneGroups(YPD_magic_t2__SlowRearranged,
                                                           'list', YPD_magic_t2_myClusteredGroups)
stat_YPD_magic_t2__SlowRearranged.to_csv('./ysc_Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged.csv')

###################################################################################################################

# 2. Do KNN for slow cells based on average values for gene set clusters
stat_YPD_magic_t2__SlowRearranged = pd.read_csv('./ysc_Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged.csv',
                                                index_col=0)
YPD_magic_t2__SlowRearranged = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2__SlowRearranged.csv', index_col=0)
phate_op = phate.PHATE()
dataPhate = phate_op.fit_transform(YPD_magic_t2__SlowRearranged)

df = getKmeansClusteringOnSelectedFeatures(stat_YPD_magic_t2__SlowRearranged,
                                                                    5, 'MAGIC_allGenes.png', 'MAGIC: all genes',
                                                                    dataPhate)
df.to_csv('./ysc_Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged_W_CellClusters.csv')


# 2A. plot dataPhate on all cells (fast included) with defined above clustering
YPD_magic_t2 = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2.csv', index_col=0)
YPD_magic_t2.drop(columns=['Cells', 'Genotype'], inplace=True)
phate_op = phate.PHATE()
dataPhate = phate_op.fit_transform(YPD_magic_t2)
df = pd.read_csv('./ysc_Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged_W_CellClusters.csv', index_col=0)
YPD_magic_t2 = YPD_magic_t2.merge(df, how='left', left_index=True, right_index=True)
YPD_magic_t2.fillna(1000, inplace=True)

viridis = cm.get_cmap('viridis', 6)
scprep.plot.scatter2d(dataPhate, c=YPD_magic_t2['clusterPredictions'], cmap=viridis)
plt.tight_layout()
plt.show()
plt.savefig('./ysc_Gresham2019/figures/phateAllCells_w_slowClusters.png')

###################################################################################################################

# 3. Jitter plots: for each cell cluster defined previously as a bin -- average across. Colors by clusters
YPD_magic_t2 = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2.csv', index_col=0)
totalStat = YPD_magic_t2.describe()
#
df = pd.read_csv('./ysc_Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged_W_CellClusters.csv', index_col=0)
YPD_magic_t2__SlowRearranged = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2__SlowRearranged.csv',
                                           index_col=0)
dfExtended = pd.merge(df, YPD_magic_t2__SlowRearranged, how='left', left_index=True, right_index=True)
dfNormalized = getLogFcExpression(dfExtended, YPD_magic_t2__SlowRearranged.columns, totalStat, '50%')

dfSlowMarkers = pd.read_csv('./ysc_Gresham2019/ExternalData/high TMRE profile.csv', sep='\t')
dfSlowMarkers = addORForGene(dfSlowMarkers)
slowMarkersDict = {}
unqGroups = list(set(dfSlowMarkers['Description']))
for i in range(len(unqGroups)):
    slowMarkersDict.update({unqGroups[i]: dfSlowMarkers[dfSlowMarkers['Description'] == unqGroups[i]]['ORF'].tolist()})

statDfNormalized = getValuesForGeneGroups(dfNormalized, 'dict', slowMarkersDict)

# get figure for Dhar VS vanDijk
externalDataAnalysis()

# get jitter
for key, value in slowMarkersDict.items():
    values = statDfNormalized[key]
    scprep.plot.jitter(statDfNormalized['clusterPredictions'], values, c=values, title=key, cmap=['red', 'blue'],
                       xlabel='cluster',
                       fontsize=14, colorbar=False, xticks=False, s=12, means_s=40, means_c='grey', figsize=(5, 5))
    plt.plot([0, 5], [1, 1], 'k--')
    plt.tight_layout()
    plt.show()
    plt.savefig('./ysc_Gresham2019/figures/jitter_SlowMarkers__CellClusters/' + key + '.png')

###################################################################################################################

# 4. in phate space plot cells but color by values in slow markers
YPD_magic_t2__SlowRearranged = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2__SlowRearranged.csv',
                                                index_col=0)
phate_op = phate.PHATE()
dataPhate = phate_op.fit_transform(YPD_magic_t2__SlowRearranged)
unqGroups = list(set(dfSlowMarkers['Description']))
fig = plt.figure(figsize=(15, 7))
i = 0
for i1 in range(len(unqGroups)):
    i += 1
    ax = plt.subplot(2, 4, i)
    plt.scatter(dataPhate[:, 0], dataPhate[:, 1], c=statDfNormalized[unqGroups[i1]], s=5, cmap='Spectral',
                norm=MidpointNormalize(midpoint=1, vmin=0.9, vmax=1.1))
    ax.set_title(unqGroups[i1], fontsize=10)
plt.colorbar()
plt.savefig('./ysc_Gresham2019/figures/scatter_PHATE_colored_SlowMarker.png')

###################################################################################################################

# 5. pairwise scatter plots between each pair of slow markers
unqGroups = list(set(dfSlowMarkers['Description']))
unqGroups.remove('Gasch RP')
viridis = cm.get_cmap('viridis', 6)
fig = plt.figure(figsize=(15, 7))
numGroupsRepresentative = len(unqGroups)
i = 0
for i1 in range(len(unqGroups)):
    for i2 in range(len(unqGroups) - 1):
        if i1 > i2:
            i += 1
            ax = plt.subplot(4, 6, i)
            X1 = statDfNormalized[unqGroups[i1]]
            X2 = statDfNormalized[unqGroups[i2]]

            plt.scatter(X1, X2, c=df['clusterPredictions'], s=5, cmap=viridis)
            plt.plot([0.9, 1.1], [1, 1], 'k--')
            plt.plot([1, 1], [0.9, 1.1], 'k--')
            ax.set_xlabel(unqGroups[i1], fontsize=8)
            ax.set_ylabel(unqGroups[i2], fontsize=8)

plt.tight_layout()
plt.savefig('./ysc_Gresham2019/figures/scatters_pairwise_SlowMarkers_w_CellClusters.png')






























