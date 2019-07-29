import scprep
import pandas as pd
import numpy as np
import magic
import phate
import math
import statsmodels.api as sm
from sklearn import linear_model


##################################################################################################################

# 1. Load a dataFrame with scRNA-seq for YPD data.
# !!!: unzip data first in your local directory

from ysc_Gresham2019.functions.preprocessGresham2019__ysc import preprocessGresham2019__ysc
YPD = preprocessGresham2019__ysc(['~/Develop/PycharmProjects/ysc/ysc_Gresham2019/ExternalData/Gresham2019/GSM3564448_YPD-fastqTomat0-Counts.tsv'],
                                 'WT(ho)')
YPD = YPD[0]
# save to .csv
YPD.to_csv('./ysc_Gresham2019/Data/YPD.csv')

##################################################################################################################

# 2. Make a YPD-MAGIC dataset

YPD = pd.read_csv('./ysc_Gresham2019/Data/YPD.csv')
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
YPD_magic_t2.to_csv('./ysc_Gresham2019/Data/YPD_magic_t2.csv')

##################################################################################################################

# 3. Make a DS containing log2FC from Riddhi's (high TMRE) and David's (slow cells):

# read Dhar's data
Dhar19 = pd.read_csv('./ysc_Gresham2019/ExternalData/Dhar2019/Dhar2019.csv', sep='\t')
Dhar19['Dhar_log_FC'] = np.log2(pd.to_numeric(Dhar19['HIGH']) / pd.to_numeric(Dhar19['LOW']) + 1)
Dhar19 = Dhar19[['Gene', 'Dhar_log_FC']]
Dhar19.replace([np.inf, -np.inf], np.nan, inplace=True)
Dhar19.dropna(inplace=True)
# read vanDijk's data
vanDijk15 = pd.read_csv('./ysc_Gresham2019/ExternalData/vanDijk2015/vanDijk2015.csv', sep='\t')
vanDijk15['vanDijk_log_FC'] = np.log2(pd.to_numeric(vanDijk15['Slow']) / pd.to_numeric(vanDijk15['Fast']) + 1)
vanDijk15 = vanDijk15[['ORF', 'Gene', 'vanDijk_log_FC']]
vanDijk15.replace([np.inf, -np.inf], np.nan, inplace=True)
vanDijk15.dropna(inplace=True)
# combine together and save
slowStat = pd.merge(Dhar19, vanDijk15, how='inner', on='Gene')
slowStat.to_csv('./ysc_Gresham2019/Data/slowStat.csv')

##################################################################################################################

# MOVE HERE A SCRIPT TO GENERATE stat_YPD_magic_t2__SlowRearranged_W_CellClusters

##################################################################################################################

# 4. Add to MAGIC YPD dataset: PHATE coordinates and clustering (based on df)

YPD_magic_t2 = pd.read_csv('./ysc_Gresham2019/Data/YPD_magic_t2.csv', index_col=0)
YPD_magic_t2.drop(columns=['Cells', 'Genotype'], inplace=True)
phate_op = phate.PHATE()
dataPhate = phate_op.fit_transform(YPD_magic_t2)
df = pd.read_csv('./ysc_Gresham2019/Data/stat_YPD_magic_t2__SlowRearranged_W_CellClusters.csv', index_col=0)
YPD_magic_t2 = YPD_magic_t2.merge(df, how='left', left_index=True, right_index=True)
YPD_magic_t2.fillna(1000, inplace=True)
YPD_magic_t2['Phate 1'] = dataPhate[:, 0]
YPD_magic_t2['Phate 2'] = dataPhate[:, 1]
YPD_magic_t2.to_csv('./ysc_Gresham2019/Data/YPD_magic_t2_W_Phate_and_cellClustering.csv')





