def getValuesForGeneGroups(df, inputType, GroupData):
    '''
    :param df: dataset with scRNAseq measurements (can be imputed dataset or original)
    :param inputType: or 'file' (for manually collected set of slow markers), or 'list' (for generated by clustering)
    :param GroupData: if inputType == 'file' - filepath; if inputType == 'list' - list
    :return: dataset statGeneGroup with columns as mean across each GeneGroup (column GeneGroup)
    '''

    GeneGroups = []
    ORFs = []
    if inputType == 'file':
        for line in open(GroupData):
            line = line.rstrip("\n")
            values = line.split("\t")
            currentORFs = []
            currentGeneGroup = values[0]
            i = 1
            while i <= len(values) - 1 and values[i] != '':
                currentORFs.append(values[i])
                i += 1
            GeneGroups.append(currentGeneGroup)
            ORFs.append(currentORFs)
    if inputType == 'list':
        for i in range(len(GroupData)):
            GeneGroups.append('Gene cluster #' + str(i))
            ORFs.append(GroupData[i])

    # make an inside function that for one cell and all GeneGroups returns mean across all ORFs in this group
    def getStatForGeneGroups(row, GeneGroups, ORFs):
        for i in range(len(GeneGroups)):
            data = row[ORFs[i]]
            row[GeneGroups[i]] = data.mean()
        return row

    statGeneGroup = df.apply(getStatForGeneGroups, args=(GeneGroups, ORFs), axis=1)

    # add which columns to save
    potentialColumnsToSave = ['Cells', 'Genotype', 'GR_rank', 'GrowthRate']
    columnsToSave = []
    for i in potentialColumnsToSave:
        if i in statGeneGroup.columns:
            columnsToSave.append(i)
    for i in range(len(GeneGroups)):
        columnsToSave.append(GeneGroups[i])

    statGeneGroup = statGeneGroup[columnsToSave[:]]

    return statGeneGroup


