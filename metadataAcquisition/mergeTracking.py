import os
import pandas as pd

def mergeTracking(dataframe, filepath):

    os.chdir(filepath)

    trackResults = pd.read_csv('timeCourseResults.txt', delimiter="\t")
    trackCellsDF = trackResults[['cell.index', 'track.index', 'cell.frame']]
    trackCellsDF.columns = ['CellNumber', 'TrackIndex', 'TimeFrame']

    trackedCellsIdentified = pd.merge(dataframe, trackCellsDF,  how='left', on=['CellNumber', 'TimeFrame'])
    trackedCellsIdentified['TimeFrame'] = trackedCellsIdentified['TimeFrame'].astype(int)

    return(trackedCellsIdentified)