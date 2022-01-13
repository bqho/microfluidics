#%%
import numpy as np
import pandas as pd


def calc_Myo1FP_Intensity(dataframe):

    trackedCells = np.unique(dataframe['TrackIndex'])
    trackedCells = trackedCells[trackedCells > 0]

    validTrackedCells = []
    Myo1ID = []
    foldChangeKO = []
    foldChangeKA = []
    KOdirection = []
    KAdirection = []

    for i in trackedCells:
        trackSubset = dataframe[dataframe['TrackIndex'] == i]

        # only accept cells that have been tracked in the experiment for
        # more than 15 frames (105 minutes)
        if len(trackSubset) > 30:
            validTrackedCells.append(i)
        else:
            continue

        mKOmax = trackSubset.iloc[trackSubset['99PercentileKO_normalized'].argmax(),:]
        mKOmin = trackSubset.iloc[trackSubset['99PercentileKO_normalized'].argmin(),:]

        if mKOmax['TimeFrame'] > mKOmin['TimeFrame']:
            KOdirection.append('Increase')
        else:
            KOdirection.append('Decrease')

        mKOsignalChange = mKOmax['99PercentileKO_normalized']/mKOmin['99PercentileKO_normalized']



        mKAmax = trackSubset.iloc[trackSubset['99PercentileKA_normalized'].argmax(),:]
        mKAmin = trackSubset.iloc[trackSubset['99PercentileKA_normalized'].argmin(),:]

        if mKAmax['TimeFrame'] > mKAmin['TimeFrame']:
            KAdirection.append('Increase')
        else:
            KAdirection.append('Decrease')

        mKAsignalChange = mKAmax['99PercentileKA_normalized']/mKAmin['99PercentileKA_normalized']

        # Now determine whether the given tracked cells is mKO of mKa Myo1
        if mKOsignalChange > mKAsignalChange:
            Myo1ID.append('Myo1_mKO')
            foldChangeKO.append(mKOsignalChange)
            foldChangeKA.append(mKAsignalChange)

        elif mKOsignalChange < mKAsignalChange:
            Myo1ID.append('Myo1_mKa')
            foldChangeKO.append(mKOsignalChange)
            foldChangeKA.append(mKAsignalChange)

        else:
            Myo1ID.append('Cannot Be Determined')


    Myo1Identity = {
        'TrackID' : validTrackedCells,
        'Myo1Identity' : Myo1ID,
        'KOFoldChange' : foldChangeKO,
        'KOdirection' : KOdirection,
        'KAFoldChange' : foldChangeKA,
        'KAdirection' : KAdirection
    }

    Myo1Identity = pd.DataFrame(Myo1Identity)

    return(Myo1Identity)



#%%
# We can really only trust the Myo1 marker reporter
# to determine which strains are which as long as there is a
# 2-fold difference between the fold-change of the KO or KA
# fluorescence between one another. Otherwise, the signals are
# too similar, and a strain may be incorrectly labeled as mKate2 or mKOk 

def myo1ID(dataframe):

    myo1Dataframe = calc_Myo1FP_Intensity(dataframe)

    TrackID = np.unique(myo1Dataframe['TrackID'])
    ID = []
    KO_KA_change = []

    for i in TrackID:
        trackIDSubset = myo1Dataframe[myo1Dataframe['TrackID']==i]
        foldChange = np.max(trackIDSubset.iloc[0,[2,4]])/np.min(trackIDSubset.iloc[0,[2,4]])
        KO_KA_change.append(foldChange)
        ID.append(i)

    Myo1Final = {
        'TrackID' : ID,
        'FoldChange' : KO_KA_change
    }

    Myo1Final = pd.DataFrame(Myo1Final)
    Myo1Final = Myo1Final[Myo1Final['FoldChange'] > 2]

    # And then merge with the original dataframe
    tracksKeep = np.unique(Myo1Final['TrackID'])
    Myo1Identity = myo1Dataframe[myo1Dataframe['TrackID'].isin(tracksKeep)]

    finalDataFrame = pd.merge(dataframe, Myo1Identity, left_on='TrackIndex', right_on='TrackID')

    return(finalDataFrame)




#
# Finally, add the division time as part of the analysis. For each
# tracked cell, determine what the doubling time calculated

# I think it's only reasonable to select only the cells that were there since
# the start of the experiment, which will make the sine curve fit most
# accurate

# myo1DataWithTracks = divTime.mergeMyo1TrackID(myo1Data, finalQuantification)

# divisiontime = []
# track = []
# count = 0

# # Subset for the indicate fluorescence channel
# myo1Fluo = 'KO'
# fluoSubset = myo1Data[myo1Data['Myo1ID']==myo1Fluo]


# for i in np.unique(fluoSubset['TrackIndex']):
#     count = count + 1
#     print(count)

#     # Subset for trackIDs that are present before MMS treatment
#     trackSubset = fluoSubset[fluoSubset['TrackIndex'] == i]
#     trackSubset = trackSubset[trackSubset['TimeFrame'] < 45]

#     # Subset for trackIDs that are present for the majority of
#     # the untreated acclimitization period
#     if len(trackSubset) > 30:

#         try:
#         # Calculate the division time
#             time = divTime.calcDivisionTime(myo1Data=myo1Data,
#                                             myo1Fluo=myo1Fluo,
#                                             trackedCells=finalQuantification,
#                                             trackIndex=i)
#             divisiontime.append(time)
#             track.append(i)

#         except:
#             pass

#     else:
#         continue

# divtime = {
#     'TrackID' : track,
#     'DivisionTime' : divisiontime
# }

# divtime = pd.DataFrame(divtime)
