def getDiffExpressedGenes(df, column, quantileThreshold):
    '''
    :param df:
    :param column:
    :param quantile:
    :return:
    '''
    N = round(quantileThreshold * len(df))

    result = df.sort_values([column], ascending=True)
    genesAscending = result.iloc[:N]['Gene'].tolist()

    result = df.sort_values([column], ascending=False)
    genesDescending = result.iloc[:N]['Gene'].tolist()

    genes = genesAscending + genesDescending
    return genes
