# ________________________________________________________________________________
import os
import pandas as pd
import glob



# ________________________________________________________________________________
# First generate the dataframe / directory for all images of interest

def imageDirectory(x):
    """[Loop through your current directory, acquire path and names of all TIF images

    The function places the path and name of files in a dataframe, where
    each row corresponds to a particular position, time frame. This dataframe is
    used as a reference for downstream analyses, providing a reference to where
    images reside, and which images to open for specific conditions]

    Args:
        x ([str]): [full path to top-level folder for image files]
    """

    # Changes to the specific image directory
    os.chdir(x)

    # images of interest in this analyses are GFP, Myo1-mKOk, and Myo1-mKate2
    imgFilesGFP = []
    imgFilesmKO = []
    imgFilesmKa = []

    # Finds all the images in the folder that are the fluorescence
    for root, dirs, files in os.walk(os.getcwd()):
        for i in files:
            if i.find('GFP') > -1:
                imgFilesGFP.append({'ImageFileName': i})
            if i.find('mKa') > -1:
                imgFilesmKa.append({'ImageFileName': i})        
            if i.find('mKO') > -1:
                imgFilesmKO.append({'ImageFileName': i})

    imgFilesGFP = pd.DataFrame(imgFilesGFP) 
    imgFilesmKO = pd.DataFrame(imgFilesmKO) 
    imgFilesmKa = pd.DataFrame(imgFilesmKa) 

    # The following define grep functions specific for the fluorescence image
    # names, to acquire the position/identity of the image and the time point
    def grepPosition(x):
        return(x[x.find('position')+len('position'):x.find('_time')])

    def grepTime(x):
        return(x[x.find('time')+len('time'):x.find('.tif')])

    # Apply the above metadata functions to the image file names
    def applyMetadata(df, fluorophore):
        df['Position'] = pd.Series(df.iloc[:,0]).apply(grepPosition)
        df['Time'] = pd.Series(df.iloc[:,0]).apply(grepTime)
        df.columns = [fluorophore, 'ChipPosition', 'TimeFrame']
        return(df)

    imgFilesGFP = applyMetadata(imgFilesGFP, fluorophore = 'GFP')
    imgFilesmKO = applyMetadata(imgFilesmKO, fluorophore = 'mKO')
    imgFilesmKa = applyMetadata(imgFilesmKa, fluorophore = 'mKA')

    finalImgDF = imgFilesGFP.merge(imgFilesmKO,on=['ChipPosition', 'TimeFrame']).merge(imgFilesmKa,on=['ChipPosition', 'TimeFrame'])
    finalImgDF['TimeFrame'] = finalImgDF['TimeFrame'].astype(int)
    finalImgDF = finalImgDF.sort_values('TimeFrame')

    return(finalImgDF)



# ________________________________________________________________________________
# First generate the dataframe / directory for all images of interest

def segmentationDirectory(x):
    """[Loops through the indicated directory to find segmentation masks
    that were outputted by CellX Segmentation and Tracking Software]

    Args:
        x ([str]): [full path to the segmentation files]
    """

    # Change to the directory with the segmentation images
    os.chdir(x)

    # Create list of segmentation files
    segFiles = []

    for i in glob.glob('mask*.mat'):
        segFiles.append({'SegmentationFile': i})

    segFiles = pd.DataFrame(segFiles)

    # For some reason, the time frame labels are outputted as
    # 5 digits, whereas the initial image files are 4 digits. We
    # will have to correct that by removing a leading 0
    def grepSegmentationTime(x):
        return(x[x.find('mask_0')+len('mask_0'):x.find('.mat')])

    segFiles['TimeFrame'] = pd.Series(segFiles.iloc[:,0]).apply(grepSegmentationTime)
    segFiles['TimeFrame'] = segFiles['TimeFrame'].astype(int)
    segFiles = segFiles.sort_values('TimeFrame')

    return(segFiles)



# ________________________________________________________________________________
# We need to generate a DF which indicates which trackID corresponds to
# which cell number in the segmented image

def trackCells(x):
    """[Provides a track ID for every cell in all timepoints, where
    each track ID represents that single cell followed over the course
    of the time-lapse imaging experiment]

    Args:
        x ([str]): [full path to the location of the timeCourseResults file
        outputted from CellX Segmentation Software]
    """
    data = pd.read_csv('timeCourseResults.txt', delimiter="\t")
    trackCellsDF = data[['cell.index', 'track.index', 'cell.frame']]
    trackCellsDF.columns = ['CellNumber', 'TrackIndex', 'TimeFrame']

    return(trackCellsDF)