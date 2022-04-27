"""
Example: Lessley Parabolas

This example is based on the parabolas used as a test case by Lessley et al. 
(2005). The effect of signal registration is explored here.

Dataset Citation:
   Lessley, D., Crandall, J., Shaw, G., Kent, R., & Funk, J. (2004). "A
      normalization technique for developing corridors from individual 
      subject responses." SAE Technical Papers.
      https://doi.org/10.4271/2004-01-0288

Copyright (c) 2022 Devon C. Hartlen

"""

from arcgen import arcgen
import numpy as np
import matplotlib.pyplot as plt
import os

if not os.path.exists('outputs'):
    os.makedirs('outputs')


# Case 1: Enabling signal normalization
# Load data from .csv and convert to list of np.arrays
inputPath = 'LessleyParabolas/LessleyParabolas.csv'
rawCsv = np.genfromtxt(inputPath, delimiter=',', encoding=None)
nSignals = rawCsv.shape[1]/2
inputSignals = []
for iSig in range(int(nSignals)):
    temp = rawCsv[:,(2*iSig):(2*iSig+2)]
    index = np.logical_not(np.isnan(temp[:,0]))
    inputSignals.append(temp[index,:])

charAvg, innCorr, outCorr, processed, debug = arcgen(
    inputSignals, 
    NormalizeSignals='on',
    nResamplePoints=100, 
    CorridorRes=100)

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Lessley Parabolas with Signal Normalization')
fig.savefig('outputs/Corridors - Lessley Parabolas with Normalization.png')

# Case 2: Disabling signal normalization. This is strongly discouraged
# Load data from .csv and convert to list of np.arrays
inputPath = 'LessleyParabolas/LessleyParabolas.csv'
rawCsv = np.genfromtxt(inputPath, delimiter=',', encoding=None)
nSignals = rawCsv.shape[1]/2
inputSignals = []
for iSig in range(int(nSignals)):
    temp = rawCsv[:,(2*iSig):(2*iSig+2)]
    index = np.logical_not(np.isnan(temp[:,0]))
    inputSignals.append(temp[index,:])

charAvg, innCorr, outCorr, processed, debug = arcgen(
    inputSignals, 
    NormalizeSignals='off',
    nResamplePoints=100, 
    CorridorRes=100)

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Lessley Parabolas without Signal Normalization')
fig.savefig('outputs/Corridors - Lessley Parabolas without Normalization.png')