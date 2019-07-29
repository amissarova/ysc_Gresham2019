def externalDataAnalysis():
    #
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    slowStat = pd.read_csv('./ysc_Gresham2019/Data/slowStat.csv', index_col=0)
    dfSlowMarkers = pd.read_csv('./ysc_Gresham2019/ExternalData/high TMRE profile.csv', sep='\t')
    dfSlowMarkers = pd.merge(dfSlowMarkers, slowStat, how='left', on='Gene')

    #
    fig = plt.figure(figsize=(14, 6))
    myColors = ['orangered', 'turquoise', 'gold', 'peru', 'Green',
                'orchid', 'navy', 'deepskyblue', 'maroon']

    myLabels = list(set(dfSlowMarkers['Description']))
    color_dict = {}
    for i in range(len(myLabels)):
        color_dict.update({myLabels[i] : myColors[i]})
    plt.scatter(slowStat['Dhar_log_FC'], slowStat['vanDijk_log_FC'], s=5, c='silver')
    plt.plot([0, 4], [1, 1], 'k--')
    plt.plot([1, 1], [0, 4], 'k--')
    plt.scatter(dfSlowMarkers['Dhar_log_FC'], dfSlowMarkers['vanDijk_log_FC'], s=30,
                color=[color_dict[i] for i in dfSlowMarkers['Description']], marker='h')
    plt.xlim([0, 4])
    plt.ylim([0, 4])

    markers = [plt.Line2D([0, 0], [0, 0], color=color, marker='o', linestyle='') for color in color_dict.values()]
    plt.legend(markers, color_dict.keys(), numpoints=1)
    plt.xlabel('Dhar2015')
    plt.ylabel('vanDijk2019')
    plt.show()
    plt.savefig('./ysc_Gresham2019/figures/Dhar_VS_vanDijk.png')

    regulationThreshDown = np.log2(1.6666)
    regulationThreshUp = np.log2(2.5)

    # get groups for Dhar
    markersDhar = ['Respiration', 'DNA damage', 'Iron deficiency DOWN', 'Iron deficiency UP', 'MSN2 target']
    markersDharType = ['Down', 'Up', 'Down', 'Up', 'Down']
    genesDhar = []
    for i in range(len(markersDhar)):
        currentGenes = dfSlowMarkers[dfSlowMarkers['Description'] == markersDhar[i]]
        if markersDharType[i] == 'Down':
            currentGenes = currentGenes[currentGenes['Dhar_log_FC'] < regulationThreshDown]['Gene'].tolist()
        elif markersDharType[i] == 'Up':
            currentGenes = currentGenes[currentGenes['Dhar_log_FC'] > regulationThreshUp]['Gene'].tolist()
        for j in currentGenes:
            genesDhar.append(j)

    markersVanDijk = ['DNA damage', 'Iron deficiency DOWN', 'Iron deficiency UP', 'MSN2 target']
    markersVanDijkType = ['Up', 'Up', 'Down', 'Up']
    genesVanDijk = []
    for i in range(len(markersVanDijk)):
        currentGenes = dfSlowMarkers[dfSlowMarkers['Description'] == markersVanDijk[i]]
        if markersVanDijkType[i] == 'Down':
            currentGenes = currentGenes[currentGenes['vanDijk_log_FC'] < regulationThreshDown]['Gene'].tolist()
        elif markersVanDijkType[i] == 'Up':
            currentGenes = currentGenes[currentGenes['vanDijk_log_FC'] > regulationThreshUp]['Gene'].tolist()
        for j in currentGenes:
            genesVanDijk.append(j)

    return genesDhar, genesVanDijk






