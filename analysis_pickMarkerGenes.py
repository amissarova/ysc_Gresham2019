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

# 1. Jitter plots: binned by cell cluster, show how different gene set clustering behaves
YPD_magic_t2 = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2.csv', index_col=0)
[YPD_magic_t2__SlowRearranged, YPD_magic_t2_corrMtrxRearranged, YPD_magic_t2_myClusteredGroups] = \
	clusteredCorrMatrix(YPD_magic_t2, 'corr_AllGenes_MAGIC.png', 'MAGIC, all genes', threshGrowthRate=0.2)

df = pd.read_csv('./ysc_Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged_W_CellClusters.csv', index_col=0)
stat_YPD_magic_t2 = getValuesForGeneGroups(YPD_magic_t2, 'list', YPD_magic_t2_myClusteredGroups)
stat_YPD_magic_t2.drop({'Cells', 'Genotype'}, axis=1, inplace=True)
stat_YPD_magic_t2['clusterPredictions'] = 'all'
dfExtended = pd.concat([df, stat_YPD_magic_t2], ignore_index=True)

for i in range(len(dfExtended.columns)-1):
	values = dfExtended[dfExtended.columns[i]]
	scprep.plot.jitter(dfExtended['clusterPredictions'], values, c=values, title=i, cmap=['red', 'blue'],
						xlabel='cluster', fontsize=14, colorbar=False, xticks=False, s=12, means_s=40, means_c='grey',
						figsize=(5, 5))
	plt.tight_layout()
	plt.show()
	plt.savefig('./ysc_Gresham2019/figures/jitter_GeneCorrClusters__CellClusters/geneClust_' + str(i))

############################################################################################################

# 2. plot mean for all cells and mean for slow cells (one dot - one gene)
YPD_magic_t2 = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2.csv', index_col=0)
[YPD_magic_t2__SlowRearranged, YPD_magic_t2_corrMtrxRearranged, YPD_magic_t2_myClusteredGroups] = \
	clusteredCorrMatrix(YPD_magic_t2, 'corr_AllGenes_MAGIC.png', 'MAGIC, all genes', threshGrowthRate=0.2)

totalStat = YPD_magic_t2.describe()
totalStat_Slow = YPD_magic_t2__SlowRearranged.describe()

clusterNums = [2, 7, 11, 16]
fig = plt.figure(figsize=(15, 7))
i = 0
for i1 in range(len(clusterNums)):
	i += 1
	ax = plt.subplot(2, 2, i)
	plt.scatter(totalStat.loc['mean'][YPD_magic_t2_myClusteredGroups[clusterNums[i1]]],
				totalStat_Slow.loc['mean'][YPD_magic_t2_myClusteredGroups[clusterNums[i1]]], s=25, color='Green')
	ax.set_title('Gene cluster #' + str(clusterNums[i1]), fontsize=12)
	ax.set_xlabel('All cells')
	ax.set_ylabel('Slow cells')
plt.tight_layout()
plt.show()
plt.savefig('./ysc_Gresham2019/figures/averageExpressionForCellClusters.png')

############################################################################################################

# 3. Choose 3 the most highest expressed genes from each cluster:
for i in range(len(clusterNums)):
	sortedStat = totalStat[YPD_magic_t2_myClusteredGroups[clusterNums[i]]].sort_values(axis=1, by='mean', ascending=False)
	with open('./ysc_Gresham2019/Data/sorted_cluster_' + str(clusterNums[i]) + '.csv', 'w') as myfile:
		wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
		wr.writerow(sortedStat.columns.tolist())

	print('Gene cluster # ' + str(clusterNums[i]))
	print(sortedStat.columns[0:5])

############################################################################################################

#4. Plot scatters for each pair from cluster #2 against cluster #7
YPD_magic_t2 = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2.csv', index_col=0)
totalStat = YPD_magic_t2.describe()
#
df = pd.read_csv('./ysc_Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged_W_CellClusters.csv', index_col=0)
YPD_magic_t2__SlowRearranged = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2__SlowRearranged.csv',
                                           index_col=0)
dfExtended = pd.merge(df, YPD_magic_t2__SlowRearranged, how='left', left_index=True, right_index=True)
dfNormalized = getLogFcExpression(dfExtended, YPD_magic_t2__SlowRearranged.columns, totalStat, '50%')

ORFs_2 = ['YGR192C', 'YHR174W', 'YLR044C', 'YKL060C', 'YBL092W']
ORFs_7 = ['YLR110C', 'YDR524C-B', 'YGL123W', 'YDR025W', 'YKL006W']
fig = plt.figure(figsize=(15, 8))
viridis = cm.get_cmap('viridis', 6)
i = 0
for i1 in range(5):
	for i2 in range(5):
		i += 1
		ax = plt.subplot(5, 5, i)
		plt.scatter(dfNormalized[ORFs_2[i1]], dfNormalized[ORFs_7[i2]],
						c=df['clusterPredictions'], s=5, cmap=viridis)
		ax.set_xlabel(ORFs_2[i1], fontsize=8)
		ax.set_ylabel(ORFs_7[i2], fontsize=8)
plt.tight_layout()
plt.show()
plt.savefig('./ysc_Gresham2019/figures/representsFromClusterGroupsToSplitCellClusters.png')

# 5. mean of all cells (X) -- ratio between
YPD_magic_t2 = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2.csv', index_col=0)
totalStat = YPD_magic_t2.describe()
#
df = pd.read_csv('./ysc_Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged_W_CellClusters.csv', index_col=0)
YPD_magic_t2__SlowRearranged = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2__SlowRearranged.csv',
                                           index_col=0)
dfExtended = pd.merge(df, YPD_magic_t2__SlowRearranged, how='left', left_index=True, right_index=True)
dfNormalized = getLogFcExpression(dfExtended, YPD_magic_t2__SlowRearranged.columns, totalStat, '50%')

# for each cell cluster: make statDf
statDfNormalized = []
for i in range(5):
	currentDfNormalized = dfNormalized[dfNormalized['clusterPredictions'] == i]
	statDfNormalized.append(currentDfNormalized.describe())
fig = plt.figure(figsize=(15, 8))

i = 0
for i1 in range(5):
	for i2 in range(5):
		i += 1
		ax = plt.subplot(5, 5, i)
		plt.scatter(statDfNormalized[i1].loc['mean'][YPD_magic_t2_myClusteredGroups[2]],
					statDfNormalized[i2].loc['mean'][YPD_magic_t2_myClusteredGroups[2]],
					c=totalStat.loc['mean'][YPD_magic_t2_myClusteredGroups[2]], s=1, cmap='ocean')
