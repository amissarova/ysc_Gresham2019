def clusteredCorrMatrix(df, figname, figtitle, threshGrowthRate=1, columnsToExclude=[]):
    '''
    Calculates cross-correlation between selected features (genes or gene sets) across slow cells
    :param df: dataset for which to calculate the correlation matrix
    :param figname: title for heatmap figure
    :param figtitle: savename for heatmap figure
    :param threshGrowthRate: subset of cells by growth rate for which to calculate the correlation matrix (default=1 i.e. all cells)
    :param columnsToExclude: additional columns to exclude (if exist, default=[])
    :param threshGrowthRate: subset of cells by growth rate for which to calculate the correlation matrix (default=1 i.e. all cells)

    :return: rearranged dataset (only for selected cells) and correlation matrix. Also prints out heatmap.
    '''
    import numpy as np
    import scipy.cluster.hierarchy as spc
    import seaborn as sns
    import matplotlib.pyplot as plt

    from Gresham2019.functions.addGrowthRate import addGrowthRate
    # add growthRate (if doesn't exist already on this dataset)
    if not 'GrowthRate' in df.columns:
        df = addGrowthRate(df)

    # get the columns on which the correlation calculation should be excluded
    # type 1: columns which don't contribute to genes or gene sets
    generalColumnsToExclude = ['Cells', 'Genotype', 'Genotype_Group', 'GR_rank', 'GrowthRate']
    generalColumnsToExclude = generalColumnsToExclude + list(set(columnsToExclude) -
                                                             set(generalColumnsToExclude))
    # type 2: gene or geneSets that have 0 for all the cells
    genesAllZeros = df.columns[(df == 0).all()].tolist()
    # exclude those columns
    allColumnsToExclude = generalColumnsToExclude + genesAllZeros
    allColumns = list(df.columns.values)
    groupsForCorr = list(set(allColumns) - set(allColumnsToExclude))


    # get subset of cells with growth rate lower than threshGrowthRate (which is quantile, between 0 and 1)
    assert ('GrowthRate' in df.columns), 'add Growth Rate to filter by growth for correlation analysis'
    dfSlow = df[df['GrowthRate'] <= df['GrowthRate'].quantile(threshGrowthRate)]
    dfSlow = dfSlow[groupsForCorr]

    # calculate correlation matrix, calculate hierarchy and rearrnage matrix based on hierarchy
    corrMtrx = dfSlow.corr().values
    d = spc.distance.pdist(corrMtrx, metric='correlation')
    L = spc.linkage(d, method='complete')

    # get the clusters list
    ind = spc.fcluster(L, 0.5 * d.max(), 'distance')
    myClusters = list(set(ind))
    myClusteredGroups = []
    for i in myClusters:
        tempCluster = []
        currentIdxBool = ind == i
        currentIdxList = [j for j in range(len(currentIdxBool)) if currentIdxBool[j] == True]
        for j in currentIdxList:
            tempCluster.append(groupsForCorr[j])
        myClusteredGroups.append(tempCluster)

    # rearrange columns for corrMtrx
    columns = [dfSlow.columns.tolist()[j] for j in list((np.argsort(ind)))]
    dfSlow_rearranged = dfSlow.reindex_axis(columns, axis=1)
    corrMtrx_rearranged = dfSlow_rearranged.corr().values

    # plot the heatmap for the correlation
    sns.set()
    fig, ax = plt.subplots(figsize=(15, 7))
    # if there are not many groups (i.e. < 50), add labels and put the gap between squares; otherwise not
    if len(corrMtrx) < 50:
        labels = []
        for i in range(len(groupsForCorr)):
            labels.append(str(i) + '. ' + dfSlow_rearranged.columns[i])
        ax = sns.heatmap(corrMtrx_rearranged, linewidths=.5, yticklabels=labels, cmap='PiYG', center=0, ax=ax)
        ax.set_yticklabels(labels, fontsize=7)
    else:
        ax = sns.heatmap(corrMtrx_rearranged, cmap='PiYG', center=0, ax=ax)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
    ax.set_title(figtitle)
    plt.tight_layout()
    plt.savefig('./Gresham2019/figures/' + figname)

    # return rearranged dataset (based on hierarchy) and correlation matrix
    return dfSlow_rearranged, corrMtrx_rearranged, myClusteredGroups
