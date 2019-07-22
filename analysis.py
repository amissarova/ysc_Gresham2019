
# import general modules
import scprep
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import phate
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import csv

# import my functions
from Gresham2019.functions.getValuesForGeneGroups import getValuesForGeneGroups
from Gresham2019.scripts.clusteredCorrMatrix import clusteredCorrMatrix
from Gresham2019.scripts.getKmeansClusteringOnSelectedFeatures import getKmeansClusteringOnSelectedFeatures

############################################################################################################

# 1. For MAGIC data: make correlation plot for all genes and get the list of gene clusters
YPD_magic_t2 = pd.read_csv('./Gresham2019/Data/YPD_magic_t2.csv', index_col=0)
[YPD_magic_t2__SlowRearranged, YPD_magic_t2_corrMtrxRearranged, YPD_magic_t2_myClusteredGroups] = \
    clusteredCorrMatrix(YPD_magic_t2, 'corr_AllGenes_MAGIC.png', 'MAGIC, all genes', threshGrowthRate=0.2)
# save rearranged dataset and clusters themselves
YPD_magic_t2__SlowRearranged.to_csv('./Gresham2019/Data/YPD_magic_t2__SlowRearranged.csv')
for i in range(len(YPD_magic_t2_myClusteredGroups)):
    with open('./Gresham2019/Data/MAGIC_clusteredGeneSets/geneSet_' + str(i) + '.csv', 'w', newline='\n') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(YPD_magic_t2_myClusteredGroups[i])
allGroupedGenes = YPD_magic_t2__SlowRearranged.columns
with open('./Gresham2019/Data/MAGIC_clusteredGeneSets/allGenes.csv', 'w', newline='\n') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(allGroupedGenes)

YPD_magic_t2__SlowRearranged = pd.read_csv('./Gresham2019/Data/YPD_magic_t2__SlowRearranged.csv', index_col=0)
# make and save stat dataset: for each gene cluster calculate the average across all the genes in the cluster
stat_YPD_magic_t2__SlowRearranged = getValuesForGeneGroups(YPD_magic_t2__SlowRearranged,
                                                           'list', YPD_magic_t2_myClusteredGroups)
stat_YPD_magic_t2__SlowRearranged.to_csv('./Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged.csv')


# 2. Do KNN for slow cells based on average values for gene set clusters

[df, variationExplained_df] = getKmeansClusteringOnSelectedFeatures(stat_YPD_magic_t2__SlowRearranged,
                                                                    6, 'MAGIC_allGeneSets.png', 'MAGIC: all gene sets')

df.to_csv('./Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged_W_CellClusters.csv')


# 3. Jitter plots: binned by cell cluster, show how different gene set clustering behaves

stat_YPD_magic_t2 = getValuesForGeneGroups(YPD_magic_t2, 'list', YPD_magic_t2_myClusteredGroups)
stat_YPD_magic_t2.drop({'Cells', 'Genotype'}, axis=1, inplace=True)
stat_YPD_magic_t2['clusterPredictions'] = 'all'
dfExtended = pd.concat([df, stat_YPD_magic_t2], ignore_index=True)

for i in range(len(dfExtended.columns)):
    values = dfExtended[dfExtended.columns[i]]
    scprep.plot.jitter(dfExtended['clusterPredictions'], values, c=values, title=i, cmap=['red', 'blue'], xlabel='cluster',
                       fontsize=14, colorbar=False, xticks=False, s=12, means_s=40, means_c='grey', figsize=(5, 5))
    plt.tight_layout()
    plt.show()
    plt.savefig('./Gresham2019/figures/jitter_GeneCorrClusters__CellClusters/geneClust_' + str(i))


# 4. Pairwise scatter plots for gene set clusters (the ones where there is a GO enrichment with GO)

groupsRepresentative = ['Gene cluster #1', 'Gene cluster #2', 'Gene cluster #4', 'Gene cluster #6',
                        'Gene cluster #7','Gene cluster #11', 'Gene cluster #16']

labels = ['#1', '#2', '#4', '#6', '#7', '#11', '#16']
fig = plt.figure(figsize=(15, 7))
numGroupsRepresentative = len(groupsRepresentative)
i = 0
for i1 in range(len(groupsRepresentative)):
    for i2 in range(len(groupsRepresentative) - 1):
        if i1 > i2:
            i += 1
            ax = plt.subplot(3, 7, i)
            X1 = df[groupsRepresentative[i1]]
            X2 = df[groupsRepresentative[i2]]

            plt.scatter(X1, X2, c=df['clusterPredictions'], s=5, cmap='tab10')
            plt.xticks([])
            plt.yticks([])
            ax.set_xlabel(labels[i1], fontsize=12)
            ax.set_ylabel(labels[i2], fontsize=12)
plt.tight_layout()
plt.savefig('./Gresham2019/figures/scatters_GeneCorrClusters_wCellClusters.png')


# 5. Individual genes against the gene clusters they correspond to (with colored cell clusters)

dfExtended = pd.merge(df, YPD_magic_t2__SlowRearranged, how='right', left_index=True, right_index=True)
myClusters = ['Gene cluster #2', 'Gene cluster #1', 'Gene cluster #7', 'Gene cluster #7',
              'Gene cluster #2', 'Gene cluster #2', 'Gene cluster #2']
myORFs = ['YML100W', 'YML032C', 'YDR263C', 'YGL013C', 'YBL005W', 'YMR037C', 'YML032C']
myGenes = ['TSL1', 'RAD52', 'DIN7', 'PDR1', 'PDR3', 'MSN2', 'MSN4']

fig = plt.figure(figsize=(15, 7))
i = 0
for i1 in range(len(myClusters)):
    i += 1
    ax = plt.subplot(2, 4, i)
    plt.scatter(dfExtended[myClusters[i1]], dfExtended[myORFs[i1]], c=dfExtended['clusterPredictions'], s=5, cmap='tab10')
    ax.set_title(myGenes[i1])
plt.savefig('./Gresham2019/figures/scatters_RepresentativeGenes.png')

# 6. Run jitter for RP genes: check if they differ by growth rate

# 7. Based on gene sets discovered from analysis of MAGIC data: do we see proper clustering for raw YPD data

















