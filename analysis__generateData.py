import scprep
import pandas as pd
import numpy as np
import magic
import statsmodels.api as sm
from sklearn import linear_model


#######################################################################
# 1. load a dataFrame with scRNA-seq for YPD data.
# save: './Gresham2019/Data/YPD.csv'
# !!! Don't forget that filename should be a list (even if contains one entry)
from Gresham2019.functions.preprocessGresham2019__ysc import preprocessGresham2019__ysc
YPD = preprocessGresham2019__ysc(['~/Develop/PycharmProjects/ysc/Gresham2019/Externaldata/Gresham2019/GSM3564448_YPD-fastqTomat0-Counts.tsv'],
                                 'WT(ho)', '~/Develop/PycharmProjects/ysc/Gresham2019/Externaldata/Gresham2019/GSM3564448_YPD_genes.tsv')
YPD = YPD[0]
# save to .csv
YPD.to_csv('./Gresham2019/Data/YPD.csv')

############################################################################

# 2. make a YPD-MAGIC dataset
YPD = pd.read_csv('./Gresham2019/Data/YPD.csv')
genes = list(YPD.columns.values)
genes.remove('Cells')
genes.remove('Genotype')
genes.remove('Genotype_Group')

# do MAGIC with t = 2
magic_operator = magic.MAGIC(t=2)
YPD_magic_t2 = magic_operator.fit_transform(YPD[genes], genes=genes)

# add 'Cells' and 'Genotype' columns from YPD
YPD_magic_t2 = pd.merge(YPD_magic_t2, YPD[['Cells', 'Genotype']], how='outer', right_index=True, left_index=True)

# save dataset
YPD_magic_t2.to_csv('./Gresham2019/Data/YPD_magic_t2.csv')








#######################################################################
#######   SUUPP
# SUPP1. Make a dataFrame where every column represents average across GOs (for GOs with no less than 25 entries).
# save: './Gresham2019/Data/GOStat_moreThan25ORFs_WT.csv'.
YPD = pd.read_csv('./Gresham2019/Data/YPD.csv')

# a) make a GO-dictionary, where GO is a key and all ORFs belonging to this GO is a value
GO = []
ORFs = []
for line in open('./Gresham2019/Externaldata/scerevisiae.GO_BP_CC_MF.ENSG_modify.gmt'):
    values = line.split("\t")
    currentORFs = []
    currentGO = values[0]
    i = 2
    while values[i] != '' and i < len(values)-1:
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
# c) apply the function to everything
GOStat = YPD.apply(getStatForGO, args=(GO, ORFs), axis=1)
GOStat = GOStat[columnsToSave[1:]]
# d) save to .csv
GOStat.to_csv('./Gresham2019/Data/GOStat_moreThan25ORFs_WT.csv')


#########
# SUPP 2. Make a stat dataset for each cell and for 3 groups from Gasch2017
# save to './Gresham2019/Data/GaschStat.csv'.

YPD = pd.read_csv('./Gresham2019/Data/YPD.csv')
# a) make an ANNO dataset: which genes belong to RP, RiBi and iESR
ANNO = pd.read_excel('~/Develop/PycharmProjects/ysc/Gresham2019/Externaldata/Gasch2017/genesToGeneFamilies.xlsx')
ANNO.replace('SGD-annotated RP', 'RP', inplace=True)
ANNO.replace('RP cluster', 'RP', inplace=True)
ANNO.replace('RiBi (originally called PAC) cluster', 'RiBi', inplace=True)
ANNO.replace('iESR cluster', 'iESR', inplace=True)

# b) function that for one cell and GeneGroup returns stat for this GeneGroup
def getStatForFamily(row, ANNO, familyName):
    df = ANNO.loc[(ANNO['Gene Group'] == familyName)]
    df = df['UID']
    rows = df.values.tolist()
    data = row[row.index.intersection(rows)]
    row[familyName] = data.values.tolist()
    row['mean_' + familyName] = np.mean(data.values.tolist())
    row['median_' + familyName] = np.median(data.values.tolist())
    row['max_' + familyName] = np.max(data.values.tolist())
    row['std_' + familyName] = np.std(data.values.tolist())

    return row[{'Cells', 'Genotype', familyName, 'mean_' + familyName, 'median_' + familyName, 'max_' + familyName, 'std_' + familyName}]

# c) a stat dataset for each cell and for 3 groups from Gasch2017
GaschStat = pd.DataFrame()
familyNames = ['RP', 'iESR', 'RiBi']
for i in familyNames:
    tempGaschStat = YPD.apply(getStatForFamily, args=(ANNO, i), axis=1)
    if GaschStat.empty:
        GaschStat = tempGaschStat
    else:
        GaschStat = pd.merge(GaschStat, tempGaschStat, how='outer', on=['Cells', 'Genotype'])
# d) save to .csv
GaschStat.to_csv('./Gresham2019/Data/GaschStat.csv')
