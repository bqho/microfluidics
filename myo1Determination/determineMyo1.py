#%%
import numpy as np
import pandas as pd

# Assign the Myo1 identity for each track index based on the 
# looking at which signal (KO or KA) had the greatest fold change
# over time

def myo1Assignment(quantifiedMyo1Images, trackedCells):
    """[this function takes the quantified pixel intensities from the
    mKO and mKate intensity images for each single cell, and determines
    whether or not the cell is likely expressing mKO or mKate2 fluorophore]

    Args:
        quantifiedMyo1Images ([dataframe]): [dataframe to the quantified mKO and mKA images]
        trackedCells ([dataframe]): [dataframe contaiing the TrackIndex, TimePoint, and
        Cell Number, outputted from CellX software]
    """
    # First merge the quanitifed myo1 intensities with tracked cells dataframe
    # to link which cell in each frame corresponds to which cell tracking
    quantifiedCells_wTracking = pd.merge(quantifiedMyo1Images, trackedCells, on =['CellNumber', 'TimeFrame'])

    # Create lists to which we will fill the single cell measured parameters
    TrackIndex = []
    Myo1ID = []
    mKAratio = []
    mKOratio = []

    count = 0

    # All the unique track indices in the experiment
    uniqueTracks = np.unique(quantifiedCells_wTracking['TrackIndex'])

    for trackID in uniqueTracks:

        # Keeping track of which tracked cell we are one
        count = count + 1
        print(str(count) + " out of " + str(len(uniqueTracks)))

        # Subset for each specified tracked cell
        track = quantifiedCells_wTracking[quantifiedCells_wTracking['TrackIndex'] == trackID]

        # Determine the maximal and minimal intensities for each
        # fluorescent channel across the time series 
        mKO = np.percentile(track['99PercentileKO'], 99)/np.percentile(track['99PercentileKO'], 1)
        mKA = np.percentile(track['99PercentileKA'], 99)/np.percentile(track['99PercentileKA'], 1)

       # Assign either the track index gets a KO or KA designation
       # based on which fold change intensity is greatest 
        if mKO > mKA:
            Myo1ID.append('KO')
        else:
            Myo1ID.append('KA')

        # Fill lists with the appropriate parameters
        TrackIndex.append(trackID)
        mKAratio.append(np.log2(mKA))
        mKOratio.append(np.log2(mKO))


    # Generate the data frame that will be used by the rest of the analysis
    myo1Assigned = {
        'TrackIndex' : TrackIndex,
        'Myo1ID' : Myo1ID,
        'log2_KA_change' : mKAratio,
        'log2_KO_change' : mKOratio,
    }

    myo1Assigned = pd.DataFrame(myo1Assigned)
    return(myo1Assigned)