"""
Example: NBDL 15g Frontal Acceleration Load - Head Acceleration

This example computes corridors for x and z linear acceleration as well as y
rotational acceleration. This example also uses three input formats: specifying
a path to read from, a list of numpy arrays, and a list of dictionaries to 
specify specimen ID

Dataset Citation:
   Ewing, C. L., & Thomas, D. J. (1972). "Human Head and Neck Response to
      Impact Acceleration." Naval Aerospace Medical Research Lab
      Pensacola Fl.

   National Highway Traffic Safety Administration. (2017). "Biomechanics
      Test Database." 
      https://www.nhtsa.gov/research-data/databases-and-software

Copyright (c) 2022 Devon C. Hartlen

"""

from arcgen import arcgen
import numpy as np
import matplotlib.pyplot as plt
import os

if not os.path.exists('outputs'):
    os.makedirs('outputs')

# Head X acceleration
# This test cases uses a path to a .csv file. Uses 4 warping points and a
# a penalty factor of 1e-2
inputPath = 'NBDL_15gFrontalDeceleration/NBDL_15gFrontal_HeadXAccel.csv'
charAvg, innCorr, outCorr, processed, debug = arcgen(
    inputPath, 
    nResamplePoints=500, 
    CorridorRes=500, 
    nWarpCtrlPts=4,
    resultsToFile=False,
    Diagnostics='off')

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('Time (ms)')
plt.ylabel('Acceleration (g)')
plt.title('X Acceleration')
fig.savefig('outputs/Corridors - X-Acceleration.png')

# Head Z acceleration
# This example uses a list of np.arrays as input to ARCGen. This is a useful
# format for when signals are preprocessed prior to computing corridors. Z-
# acceleration uses 4 control points and a penalty factor of 1e-2

# Load data from .csv and convert to list of np.arrays
inputPath = 'NBDL_15gFrontalDeceleration/NBDL_15gFrontal_HeadZAccel.csv'
rawCsv = np.genfromtxt(inputPath, delimiter=',', encoding=None)
nSignals = rawCsv.shape[1]/2
inputSignals = []
for iSig in range(int(nSignals)):
    temp = rawCsv[:,(2*iSig):(2*iSig+2)]
    index = np.logical_not(np.isnan(temp[:,0]))
    inputSignals.append(temp[index,:])

charAvg, innCorr, outCorr, processed, debug = arcgen(
    inputSignals, 
    nResamplePoints=500, 
    CorridorRes=500, 
    nWarpCtrlPts=4)

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('Time (ms)')
plt.ylabel('Acceleration (g)')
plt.title('Z Acceleration')
fig.savefig('outputs/Corridors - Z-Acceleration.png')

# Head Y Rotational Accelration
# This example uses a list of dictionaries as input. Useful when one wants to
# specify the specimen ID as well as data. Y rotational acceleration uses 5
# warping control points and a penalty factor of 1e-2

# Load data from .csv and convert to list of np.arrays
inputPath = 'NBDL_15gFrontalDeceleration/NBDL_15gFrontal_HeadYRotAccel.csv'
# Load header information from file
with open(inputPath) as myfile:
    head = [next(myfile) for x in range(1)]
head = head[0].rstrip('\n')
head = head.split(',')
specIds = []
for entry in head:
    if len(entry) > 0:
        specIds.append(entry)

# Load data itself
rawCsv = np.genfromtxt(inputPath, 
    delimiter=',',
    skip_header=1,
     encoding=None)
nSignals = rawCsv.shape[1]/2
inputSignals = []

# Construct dictionary
for iSig in range(int(nSignals)):
    temp = rawCsv[:,(2*iSig):(2*iSig+2)]
    index = np.logical_not(np.isnan(temp[:,0]))
    inputSignals.append({'data': temp, 'specId': specIds[iSig]})

charAvg, innCorr, outCorr, processed, debug = arcgen(
    inputSignals, 
    nResamplePoints=500, 
    CorridorRes=500, 
    nWarpCtrlPts=4)

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('Time (ms)')
plt.ylabel('Acceleration (rad/s^2)')
plt.title('Y Rotational Acceleration')
fig.savefig('outputs/Corridors - Y-Rot-Acceleration.png')