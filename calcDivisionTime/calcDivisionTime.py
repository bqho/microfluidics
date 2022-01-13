#%%
# Import necessary packages
import numpy as np
import pandas as pd
import os
import seaborn as sns
import scipy.optimize
import pylab as plt
import matplotlib.pyplot as plot

# Import custom package for reading metadata and other
# file I/O
import sys,inspect
current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir) 

import metadataAcquisition.data_loader as meta



# %%
# Assign the track index to each of the cells

def mergeMyo1TrackID(dataframe, trackIndex):
    trackMyo1 = pd.merge(dataframe, trackIndex, on = ['CellNumber', 'TimeFrame'])
    trackMyo1 = trackMyo1[['CellNumber', 'TrackIndex', 'ImageName', 'TimeFrame',
                           '99PercentileKO', '99PercentileKA', 'Myo1ID']]
    return(trackMyo1)



#%%
# Solution for fitting sine function obtained from
# stacked overflow, user 'unsym'from
# 'https://stackoverflow.com/questions/16716302/how-do-i-fit-a-sine-curve-to-my-data-with-pylab-and-numpy'

def fit_sin(tt, yy):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    tt = np.array(tt)
    yy = np.array(yy)
    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(np.fft.fft(yy))
    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c
    popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}



# %%
# Fit sine curve to each track index, only for
# the time points before MMS addition

def calcDivisionTime(myo1Data, myo1Fluo, trackedCells, trackIndex):
    # merge dataframe with tracking indices to identify
    # same cells tracked over time
    # myo1Tracks = mergeMyo1TrackID(myo1Data,trackedCells)

    # Subset for the specific tracked cell of interest
    myo1SingleTrack = myo1Data[myo1Data['TrackIndex']==trackIndex]

    # Subset data only for time points before MMS addition
    trackedMyo1Fluo = myo1SingleTrack[myo1SingleTrack['TimeFrame'] < 40]

    # extract relevant parameters for the sine fit
    time = trackedMyo1Fluo['TimeFrame']

    # select the appropriate Myo1 fluorescence channel, as 
    # determine by the 'determineMyo1' module in other folder
    
    if myo1Fluo == 'KO':
        myo1FluoSubset = trackedMyo1Fluo[trackedMyo1Fluo['Myo1ID']=='KO']
        signal = myo1FluoSubset['99PercentileKO']
    else:
        myo1FluoSubset = trackedMyo1Fluo[trackedMyo1Fluo['Myo1ID']=='KA']
        signal = myo1FluoSubset['99PercentileKA']

    # Fit the sine curve with function pre defined above
    res = fit_sin(time, signal)
    # print( "Amplitude=%(amp)s, Angular freq.=%(omega)s, phase=%(phase)s, offset=%(offset)s, Period.=%(period)s" % res )

    # optional: plot the curves to make sure that things make sense
    # and that the Myo1 marker is indeed oscillating
    plt.plot(time, signal, "-k", label="y", linewidth=2)
    plt.plot(time, res["fitfunc"](time), "r-", label="y fit curve", linewidth=2)
    plt.legend(loc="best")
    plt.show()

    # return the period of the sine wave, roughly corresponding
    # to the doubling time of cells
    return(res['period'])



#%%
# And just testing to make sure the function works

desktopPath = "/Users/brandonho/Desktop/"
os.chdir(desktopPath)

finalQuantification = pd.read_csv("/Users/brandonho/Desktop/finalQuantification_LSM1.csv")
myo1Data=finalQuantification[['CellNumber',
                             'TimeFrame',
                             'TrackIndex', 
                             '99PercentileKO',
                             '99PercentileKA',
                             'Myo1ID']]


calcDivisionTime(myo1Data=myo1Data,
                 myo1Fluo='KO',
                 trackedCells=finalQuantification,
                 trackIndex=565)

# Looks like it does. Now we can incorporate into the full analysis

# %%
# merge dataframe with tracking indices to identify
# same cells tracked over time
# myo1Tracks = mergeMyo1TrackID(myo1Data,trackedCells)

# Subset for the specific tracked cell of interest
myo1SingleTrack = myo1Data[myo1Data['TrackIndex']==565]

# Subset data only for time points before MMS addition
trackedMyo1Fluo = myo1SingleTrack[myo1SingleTrack['TimeFrame'] < 40]

# extract relevant parameters for the sine fit
time = trackedMyo1Fluo['TimeFrame']

# select the appropriate Myo1 fluorescence channel, as 
# determine by the 'determineMyo1' module in other folder

myo1FluoSubset = trackedMyo1Fluo[trackedMyo1Fluo['Myo1ID']=='KO']
signal = myo1FluoSubset['99PercentileKO']


# Fit the sine curve with function pre defined above
res = fit_sin(time, signal)
# print( "Amplitude=%(amp)s, Angular freq.=%(omega)s, phase=%(phase)s, offset=%(offset)s, Period.=%(period)s" % res )

# optional: plot the curves to make sure that things make sense
# and that the Myo1 marker is indeed oscillating
plt.plot(time, signal, "-k", label="y", linewidth=2)
plt.plot(time, res["fitfunc"](time), "r-", label="y fit curve", linewidth=2)
plt.legend(loc="best")

plot.savefig('trackmyo1.pdf')
# %%
