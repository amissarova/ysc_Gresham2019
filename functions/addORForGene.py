def addORForGene(df):
    '''
    :param df: 
    :return: 
    '''
    import pandas as pd
    ORF_Gene_dict = pd.read_csv('./ysc_Gresham2019/ExternalData/ORF_Gene_dict.csv', sep='\t')

    dfNew = df[:]
    if 'Gene' in dfNew.columns:
        for i, rows in dfNew.iterrows():
            index = ORF_Gene_dict[ORF_Gene_dict['Gene'] == dfNew.loc[i]['Gene']].index.tolist()
            dfNew.at[i, 'ORF'] = ORF_Gene_dict.iloc[index[0]]['ORF']

    elif 'ORF' in dfNew.columns:
        for i in range(len(dfNew)):
            index = ORF_Gene_dict[ORF_Gene_dict['ORF'] == dfNew.loc[i]['ORF']].index.tolist()
            dfNew.at[i, 'Gene'] = ORF_Gene_dict.iloc[index[0]]['Gene']

    return dfNew

