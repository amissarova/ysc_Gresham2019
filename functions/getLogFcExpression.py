def getLogFcExpression(df, columns, normDf, normMetric):
    '''
    :param df:
    :param columns:
    :param normDf:
    :param normMetric:
    :return:
    '''
    import pandas as pd
    import numpy as np
    def getLogFcExpressionForOneColumn(df, normDf, normMetric):
        currentColumn = df.name
        normFactor = normDf.loc[normMetric][currentColumn]
        newDf = df.apply(lambda x: np.log2(max(x, 0) / normFactor + 1))
        return newDf

    newDf = df[columns]
    newDf = newDf.apply(getLogFcExpressionForOneColumn, args=(normDf, normMetric), axis=0)

    if 'clusterPredictions' in df.columns:
        df = df['clusterPredictions'].to_frame()
        newDf = pd.merge(newDf, df, how='outer', left_index=True, right_index=True)

    return newDf








