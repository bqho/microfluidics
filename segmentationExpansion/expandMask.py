#%%

# Function to expand the cell segmentation masks obtained from
# CellX software, since the segmentation does not capture the 
# plamsa membrane

import pandas as pd
import numpy as np
import scipy as scipy
import os
import scipy
import scipy.io as scipyio



# %%

# segmentationDict = scipy.io.loadmat("/Users/brandonho/Documents/GitHub/brownlab-files/Python/MicrofluidicsAnalysis/sampleData/mask_00001.mat") # opens matlab .mat and stores as dictionary
# segMask = segmentationDict['segmentationMask'] # extract only segmentation mask key

# Must provide the file for which the cell expansion is
# to be performed
def expandMask(file, pixelExpansion):
    """[The function takes the cell segmentation matrix output
    from CellX Software, and expands each cell mask to increase
    the perimeter of the cell]

    Args:
        file ([string]): [the file name / path name for
        the segmentation matrix]
        pixelExpansion ([float]): [the number of pixels the user
        wants to expand the segmentation by]
    """

    # file represents the segmentation matrix that is outputted
    # by CellX software. Make a copy of the matrix in case there
    # are any issues that may override the original segmentation
    refinedSegmentation = file.copy()

    # User defines the number of pixels to expand the segmentation
    # mask by
    pixelExpansion = pixelExpansion

    # Expand the segmentation mask for all cells in the matrix
    for cellNumber in np.unique(refinedSegmentation):
        
        # Ignore the zero index
        if cellNumber == 0:
            continue
        
        print("Expanding cell #: " + str(cellNumber))

        try:
            # Identify the location / indices for where the specified cells
            # are located in the matrix
            cellIndices = np.nonzero(refinedSegmentation == cellNumber)

            # # First determine if the cell is on the border
            # if np.where(cellIndices[0]==cellNumber) < 4 :
            #     minusRow = cellIndices[0]

            # Adjust the segmentation matrix by expanding the indices
            # in the minus direction by the user specified pixel number
            minusRow = cellIndices[0] - pixelExpansion
            minusCol = cellIndices[1] - pixelExpansion

            # Relabel the new indices of that cell expansion as the current cell
            refinedSegmentation[minusRow, minusCol] = cellNumber

            # Repeat the index expansion for the current cell in the 
            # plus direction
            plusRow = cellIndices[0] + pixelExpansion
            plusCol = cellIndices[1] + pixelExpansion
            refinedSegmentation[plusRow, plusCol] = cellNumber
        except:
            continue

    return(refinedSegmentation)

# %%
