import scprep
import matplotlib as plt
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

#####################################################################################################

# This branch calculates average value of the expression for each GO and t

from Gresham2019.functions.preprocessGresham2019__ysc import preprocessGresham2019__ysc
YPD = preprocessGresham2019__ysc(['~/Develop/PycharmProjects/ysc/Gresham2019/Externaldata/Gresham2019/GSM3564448_YPD-fastqTomat0-Counts.tsv'],
                                 'WT(ho)', '~/Develop/PycharmProjects/ysc/Gresham2019/Externaldata/Gresham2019/GSM3564448_YPD_genes.tsv')
YPD = YPD[0]

# 2. Make a dataFrame where every column represents average across GOs (for GOs with no less than 25 entries).
# save: './Gresham2019/Data/GOStat_moreThan25ORFs_WT.csv'.

# a) make a GO-dictionary, where GO is a key and all ORFs belonging to this GO is a value
GO = []
ORFs = []
for line in open('./Gresham2019/Externaldata/scerevisiae.GO_BP_CC_MF.ENSG_modify.gmt'):
    values = line.split("\t")
    currentORFs = []
    currentGO = values[0]
    i = 2
    while i < len(values)-1 and values[i] != '':
        currentORFs.append(values[i])
        i += 1
    # take only GOs with >= 25 ORFs
    if len(currentORFs) >= 25 and currentGO[0:2] == 'BP':
        GO.append(currentGO)
        ORFs.append(currentORFs)

# b) function that for one cell and all GOs returns mean/median stat for all GOs
def getStatForGO(row, GO, ORFs):
    for i in range(len(GO)):
        data = row[ORFs[i]]
        row['mean_' + GO[i]] = data.mean()
    return row
columnsToSave = ['Cells', 'Genotype']
for i in range(len(GO)):
    columnsToSave.append('mean_' + GO[i])

# c) apply the function to every cell
GOStat = YPD.apply(getStatForGO, args=(GO, ORFs), axis=1)
GOStat = GOStat[columnsToSave[1:]]

#GOStat.to_csv('./Gresham2019/Data/GOStat_moreThan25ORFs_WT.csv')

#####################################################################################################

# This branch makes jitter plots: for each slow marker group (iESR, MSN2 targets, respiratory things) based on RP-rank binning

YPD_magic_t2 = pd.read_csv('./Gresham2019/Data/YPD_magic_t2.csv', index_col=0)
YPD_magic_t2 = pd.read_csv('./Gresham2019/Data/YPD.csv', index_col=0)
# 1. get values for interesting gene sets (for YPD and YPD_magic_t2)
from Gresham2019.functions.getValuesForGeneGroups import getValuesForGeneGroups
ANNO_path = './Gresham2019/Data/SlowGrowthRateMarkers.csv'
statGeneGroup = getValuesForGeneGroups(YPD, ANNO_path)
statGeneGroup_magic_t2 = getValuesForGeneGroups(YPD_magic_t2, ANNO_path)

groups = list(statGeneGroup.columns.values)
groups.remove('Cells')
groups.remove('Genotype')
groups.remove('GR_rank')

# Fig.1: Jitter plots (binned by RP) for each interesting geneSet (Raw and MAGIC)
labels = statGeneGroup['GR_rank']
labels_magic_t2 = statGeneGroup_magic_t2['GR_rank']
for i in groups:
    values = statGeneGroup[i]
    scprep.plot.jitter(labels, values, c=values, cmap=['red', 'blue'], xlabel='RP bin', title='Raw: ' + i,
                       fontsize=7, xticks=False, s=12, means_s=40, means_c='grey', figsize=(7, 7),
                       filename='./Gresham2019/figures/geneGroups_RP/' + i + '_Raw.png')
    values_magic_t2 = statGeneGroup_magic_t2[i]
    scprep.plot.jitter(labels_magic_t2, values_magic_t2, c=values, cmap=['red', 'blue'], xlabel='RP bin',
                       title='MAGIC: ' + i,
                       fontsize=7, xticks=False, s=12, means_s=40, means_c='grey', figsize=(7, 7),
                       filename='./Gresham2019/figures/geneGroups_RP/' + i + '_MAGIC.png')


############################################################################################################
# KNN clustering:
# Fig.3. Elbow method to choose number of clusters
fig = plt.figure(figsize=(10, 7))
datasets = [statGeneGroup, statGeneGroup_magic_t2]
titles = ['Raw', 'MAGIC']
threshSlowBin = 5
for i in range(len(datasets)):

    # select only slow cells (i.e. first N bins) for clustering
    idx_slowCells = datasets[i].index[datasets[i]['GR_rank'] <= threshSlowBin].tolist()
    statGeneGroup__Slow = datasets[i].iloc[idx_slowCells]

    # Fig.3. Elbow method to select num of clusters
    X = statGeneGroup__Slow[groups]
    distortions = []
    K = range(1, 15)
    for k in K:
        kmeanModel = KMeans(n_clusters=k).fit(X)
        kmeanModel.fit(X)
        distortions.append(sum(np.min(cdist(X, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / X.shape[0])

    # Plot the elbow

    ax = plt.subplot(2, 1, i + 1)
    plt.plot(K, distortions, 'kx-')
    ax.set_title(titles[i])
    ax.set_xlabel('K')
    ax.set_ylabel('Distortion')

plt.tight_layout()
plt.savefig('./Gresham2019/figures/elbowForKmeansForSlowMarkers.png')

# Fig. 4. Illustrate distortion in Raw and Magic via PHATE (for all cells and slow 20%)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 8))
list_axes = [ax1, ax2, ax3, ax4]
allGenes = list(YPD_magic_t2.columns.values)
allGenes.remove('Cells')
allGenes.remove('Genotype')
allGenes.remove('GR_rank')

datasets = [YPD, YPD_magic_t2]
titles_datasets = ['Raw: ', 'MAGIC: ']
titles_bins = ['All', 'Slow']
threshBins = [25, 5]
i = 0
for i1 in range(len(datasets)):
    phate_op = phate.PHATE()
    data_phate = phate_op.fit_transform(datasets[i1][allGenes])
    for i2 in range(len(datasets)):
        idx = datasets[i1].index[datasets[i1]['GR_rank'] <= threshBins[i2]].tolist()
        scprep.plot.scatter2d(data_phate[idx], title=titles_datasets[i1] + titles_bins[i2], ax=list_axes[i])
        i += 1
plt.tight_layout()
plt.savefig('./Gresham2019/figures/PHATE_all_and_slow.png')

# Fig 5. For MAGIC-data: cluster by 6 clusters and for each group make a jitter

numClusters = 6
threshSlowBin = 5
phate_op = phate.PHATE()
data_phate = phate_op.fit_transform(YPD_magic_t2[allGenes])
idx_slowCells = YPD_magic_t2.index[YPD_magic_t2['GR_rank'] <= threshSlowBin].tolist()

statGeneGroup_magic_t2__Slow = statGeneGroup_magic_t2[statGeneGroup_magic_t2['GR_rank'] <= threshSlowBin]
y_pred = KMeans(n_clusters=numClusters).fit_predict(statGeneGroup_magic_t2__Slow[groups])
statGeneGroup_magic_t2__Slow['clusterPredictions'] = y_pred.tolist()
statGeneGroup_magic_t2__Slow.to_csv('./Gresham2019/Data/statGeneGroup_magic_t2__Slow.csv')

scprep.plot.scatter2d(data_phate[idx_slowCells], c=y_pred, title='Slow 20%', cmap='jet')
plt.show()
plt.savefig('./Gresham2019/figures/PHATE_w_clusters.png')

# for each geneSet plot how it looks like
for i in groups:
    values = statGeneGroup_magic_t2__Slow[i]
    scprep.plot.jitter(y_pred, values, c=values, cmap=['red', 'blue'], xlabel='cluster', title='MAGIC: ' + i,
                       fontsize=7, xticks=False, s=12, means_s=40, means_c='grey', figsize=(7, 7),
                       filename='./Gresham2019/figures/Kmeans/' + i + '_MAGIC_6_clusters.png')


#
groupsRepresentative = ['Gasch:RP', 'MSN2 targets', 'PDR1,3',
                        'MF:proximal promoter sequence-specific DNA binding',
                        'Environmental DNA damage response', 'Gasch:iESR']
labels = ['RP', 'MSN2 targ.', 'PDR1,3', 'TFs', 'DNA dam. resp.', 'iESR']
fig = plt.figure(figsize=(12, 8))
numGroupsRepresentative = len(groupsRepresentative)
i = 0
for i1 in range(numGroupsRepresentative):
    for i2 in range(numGroupsRepresentative - 1):
        if i1 > i2:
            ax = plt.subplot(4, 4, i + 1)
            X1 = statGeneGroup_magic_t2__Slow[groupsRepresentative[i1]]
            X2 = statGeneGroup_magic_t2__Slow[groupsRepresentative[i2]]
            #df = X1.to_frame().join(X2.to_frame())
            #scprep.plot.scatter2d(df, c=y_pred, s=3)
            plt.scatter(X1, X2, c=y_pred, s=5, cmap='jet')
            plt.xticks([])
            plt.yticks([])
            ax.set_xlabel(labels[i1], fontsize=7)
            ax.set_ylabel(labels[i2], fontsize=7)

            i += 1
plt.tight_layout()
plt.savefig('./Gresham2019/figures/scatter_groupsRepresentativeWClusters.png')







