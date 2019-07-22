def preprocessGresham2019__ysc(filenames, strainOfInterest, geneToORFmapping):

    '''
    Reads tsv-files (one DS for one element in filenames);
    Filters out empty entries or genes with 0 coverage across all cells;
    Does library normalization and square root transfromation;
    Selects only entries for given strain of interest (e.g. only WT);
    Changes ORF names for systematic names
    Returns a list of pandas dataframes, each element -- one DS

    !!! Requires a little preprocess manual step with a filename: at the left top,
    at the corner modify the empty Excel cell and put in 'Cells'
    !!!

    !!! Don't forget that filenames should be a list (even if contains one entry)
    '''
    import scprep
    import pandas as pd

    # make a mapping between ORFs and Systematic gene names
    GeneDF = pd.Series.from_csv(geneToORFmapping, sep='\t', header=None)
    GeneDF = pd.Series.to_dict(GeneDF)


    DF = []
    for i in filenames:
        # read Gresham's data
        currentDF = scprep.io.load_csv(i, delimiter='\t')

        # skip last 5 columns: counts for barcode + ID for rows (cells))
        currentDF = currentDF.drop(labels=['KANMX', 'NATMX', 'Genotype', 'Genotype_Group', 'Replicate'], axis=1)

        # remove empty cells (rows) and genes (columns)
        currentDF = scprep.filter.filter_empty_cells(currentDF)
        #currentDF = scprep.filter.filter_empty_genes(currentDF)

        #scprep.plot.plot_library_size(YPD, cutoff=0)
        #plt.show()

        # Library size normalization and square root transformation (the latter is omitted for now)
        currentDF = scprep.normalize.library_size_normalize(currentDF)
        currentDF = scprep.transform.sqrt(currentDF)

        # use Genotype_Group and Genotype (+Replicate) to get ID vector for each cell and select only WT cells
        ID = scprep.io.load_csv(i, delimiter='\t', usecols=['Cells', 'Genotype_Group', 'Genotype'])
        ID = ID.loc[ID['Genotype_Group'] == strainOfInterest]

        # Merging YPD and ID select only WT entries from YPD
        currentDF = pd.merge(currentDF, ID, how='right', left_on='Cells', right_on='Cells')

        # changing column names from ORF to genes
        #currentDF.rename(columns=dict(zip(GeneDF.keys(), GeneDF.values())), inplace=True)

        DF.append(currentDF)

    return DF
