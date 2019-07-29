def getKmeansClusteringOnSelectedFeatures(df, numClusters, figname, figtitle, dataPhate=[]):
    '''
    :param df: dataset: rows are cells and columns are features based on which you want to do clustering
    :param numClusters: number of clusters for the cells (or manually, or using elbow method)
    :param figname: the filepath for the figure (PHATE + highlighted clusters)
    :param figtitle: the title for the figure
    :param dataPhate: 2-vector dataset for visualizing clustering (default=[] and in this case the script with calculate
    PHATE coordinates). Other possible inputs are, for example, PCA coordinates
    :return: cell clustering, and also the sign of the effect of each feature on each cluster
    '''
    import phate
    from sklearn.cluster import KMeans
    import scprep
    import matplotlib.pyplot as plt
    import pandas as pd
    import scipy.stats as stats
    import numpy as np
    from matplotlib import cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    import seaborn as sns
    # for visualization
    if dataPhate == []:
        phate_op = phate.PHATE()
        dataPhate = phate_op.fit_transform(df)

    potentialColumnsToExclude = ['Cells', 'Genotype', 'GR_rank', 'GrowthRate', 'Unnamed: 0']
    myColumns = []
    for i in df.columns:
        if i not in potentialColumnsToExclude:
            myColumns.append(i)

    # get Kmeans clustering
    y_pred = KMeans(n_clusters=numClusters).fit_predict(df[myColumns])
    df['clusterPredictions'] = y_pred.tolist()

    # plot PHATE with colored clusters
    viridis = cm.get_cmap('viridis', 6)
    scprep.plot.scatter2d(dataPhate, c=df['clusterPredictions'], title=figtitle, cmap=viridis)
    plt.tight_layout()
    plt.show()
    plt.savefig('./ysc_Gresham2019/figures/phateWithCellClusters_' + figname)

    # make a variationExplained_df: each row - cluster num; each column - one gene set.
    # value: a p-value of ttest -> values in given gene set in a given cluster VS values in given gene set in other cluster
    myClusters = list(set(y_pred))
    variationExplained_df = pd.DataFrame(myClusters, columns=['Cell cluster']).set_index('Cell cluster', drop=True)

    potentialColumnsToExclude = ['Cells', 'Genotype', 'GR_rank', 'GrowthRate', 'Unnamed: 0', 'clusterPredictions']
    myColumns = []
    for i in df.columns:
        if i not in potentialColumnsToExclude:
            myColumns.append(i)

    for i in myColumns:
        geneSetSign = []
        for j in myClusters:
            currentClusterIdx = df[df['clusterPredictions'] == j].index
            otherClustersIdx = df[df['clusterPredictions'] != j].index
            X1 = df.loc[currentClusterIdx][i]
            X2 = df.loc[otherClustersIdx][i]
            currentP = stats.ttest_ind(X1, X2).pvalue
            if currentP <= 0.05 and np.mean(X1) > np.mean(X2):
                geneSetSign.append(1)
            elif currentP <= 0.05 and np.mean(X1) < np.mean(X2):
                geneSetSign.append(-1)
            else:
                geneSetSign.append(0)
        variationExplained_df[i] = geneSetSign

    return df, variationExplained_df




















