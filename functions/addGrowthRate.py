def addGrowthRate(df):
    '''
    For each cell calculates the average across RP genes by Gasch2017 and adds the column to df
    :param df: df for which to calculate growthRate
    :return: df + growthRate
    '''
    import pandas as pd
    import numpy as np

    #growth rate will be based on RP from Gasch2017 dataset
    ANNO = pd.read_excel('~/Develop/PycharmProjects/ysc/Gresham2019/Externaldata/Gasch2017/genesToGeneFamilies.xlsx')
    ANNO.replace({'SGD-annotated RP': 'RP', 'RP cluster': 'RP', 'iESR cluster': 'iESR',
                  'RiBi (originally called PAC) cluster': 'RiBi'}, inplace=True)

    # function that for one cell and GeneGroup returns stat for this GeneGroup
    def getMeanForRP(row):
        currentGroup = ANNO.loc[(ANNO['Gene Group'] == 'RP')]
        currentGroup = currentGroup['UID']
        rows = currentGroup.values.tolist()
        data = row[row.index.intersection(rows)]

        assert (not data.empty), 'There are no RP genes in dataset to calculate growth rate. ' \
                                 'Assign growth rate to prior dataset.'
        row['GrowthRate'] = np.mean(data.values.tolist())

        return row[{'Cells', 'GrowthRate'}]

    # apply to all the cells
    growthRateStat = df.apply(getMeanForRP, axis=1)

    # merge with the main dataset
    dfWithGrowthRate = df.merge(growthRateStat, left_on='Cells', right_on='Cells')

    return dfWithGrowthRate



