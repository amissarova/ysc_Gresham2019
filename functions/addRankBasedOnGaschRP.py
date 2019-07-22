def addRankBasedOnGaschRP(df, N_bins):
    import pandas as pd
    import numpy as np

    # a) make an ANNO dataset: which genes belong to RP, RiBi and iESR
    ANNO = pd.read_excel('~/Develop/PycharmProjects/ysc/Gresham2019/Externaldata/Gasch2017/genesToGeneFamilies.xlsx')
    ANNO.replace('SGD-annotated RP', 'RP', inplace=True)
    ANNO.replace('RP cluster', 'RP', inplace=True)
    ANNO.replace('RiBi (originally called PAC) cluster', 'RiBi', inplace=True)
    ANNO.replace('iESR cluster', 'iESR', inplace=True)

    # b) function that for one cell and GeneGroup returns stat for this GeneGroup
    def getStatForFamily(row, ANNO, familyName):
        currentGroup = ANNO.loc[(ANNO['Gene Group'] == familyName)]
        currentGroup = currentGroup['UID']
        rows = currentGroup.values.tolist()
        data = row[row.index.intersection(rows)]
        row[familyName] = data.values.tolist()
        row['mean_' + familyName] = np.mean(data.values.tolist())

        return row[{'Cells', 'Genotype', familyName, 'mean_' + familyName}]

    # c) a stat dataset for each cell and for 3 groups from Gasch2017
    GaschStat = df.apply(getStatForFamily, args=(ANNO, 'RP'), axis=1)

    # d)
    labels = []
    for i in range(1, N_bins + 1):
        labels.append(i)
    GR_rank = pd.qcut(GaschStat['mean_RP'], N_bins, labels=labels)
    GaschStat = GaschStat.assign(GR_rank=GR_rank.values)

    # Add rank values to YPD dataset
    dfWithRank = df.merge(GaschStat[['Cells', 'GR_rank']], left_on='Cells', right_on='Cells')

    return dfWithRank