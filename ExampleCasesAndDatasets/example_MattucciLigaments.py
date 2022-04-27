"""
Example: Lessley Parabolas

Demonstration of ARCGen on highly variable, strictly monotonic data.

Dataset Citation:
    Mattucci, S. F. E., & Cronin, D. S. (2015). A method to characterize 
        average cervical spine ligament response based on raw data sets for 
        implementation into injury biomechanics models. Journal of the 
        Mechanical Behavior of Biomedical Materials<, 41, 251–260. 
        https://doi.org/10.1016/j.jmbbm.2014.09.023

Copyright (c) 2022 Devon C. Hartlen

"""

from arcgen import arcgen
import numpy as np
import matplotlib.pyplot as plt
import os

if not os.path.exists('outputs'):
    os.makedirs('outputs')

# Anterior Longitudinal Ligament
inputPath = 'Mattucci_CervicalLigaments/Mattucci_AnteriorLongitudinalLigament_QuasiStatic.csv'
rawCsv = np.genfromtxt(inputPath, delimiter=',', encoding=None)
nSignals = rawCsv.shape[1]/2
inputSignals = []
for iSig in range(int(nSignals)):
    temp = rawCsv[:,(2*iSig):(2*iSig+2)]
    index = np.logical_not(np.isnan(temp[:,0]))
    inputSignals.append(temp[index,:])

charAvg, innCorr, outCorr, processed, debug = arcgen(
    inputSignals, 
    nResamplePoints=200, 
    CorridorRes=100)

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('Displacement (mm)')
plt.ylabel('Force (N)')
plt.title('Anterior Longitudinal Ligament')
fig.savefig('outputs/Corridors - Anterior Longitudinal Ligament.png')

# Posterior Longitudinal Ligament
inputPath = 'Mattucci_CervicalLigaments/Mattucci_PosteriorLongitudinalLigament_QuasiStatic.csv'
rawCsv = np.genfromtxt(inputPath, delimiter=',', encoding=None)
nSignals = rawCsv.shape[1]/2
inputSignals = []
for iSig in range(int(nSignals)):
    temp = rawCsv[:,(2*iSig):(2*iSig+2)]
    index = np.logical_not(np.isnan(temp[:,0]))
    inputSignals.append(temp[index,:])

charAvg, innCorr, outCorr, processed, debug = arcgen(
    inputSignals, 
    nResamplePoints=200, 
    CorridorRes=100)

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('Displacement (mm)')
plt.ylabel('Force (N)')
plt.title('Posterior Longitudinal Ligament')
fig.savefig('outputs/Corridors - Posterior Longitudinal Ligament.png')

# Capsular Ligament
inputPath = 'Mattucci_CervicalLigaments/Mattucci_CapsularLigament_QuasiStatic.csv'
rawCsv = np.genfromtxt(inputPath, delimiter=',', encoding=None)
nSignals = rawCsv.shape[1]/2
inputSignals = []
for iSig in range(int(nSignals)):
    temp = rawCsv[:,(2*iSig):(2*iSig+2)]
    index = np.logical_not(np.isnan(temp[:,0]))
    inputSignals.append(temp[index,:])

charAvg, innCorr, outCorr, processed, debug = arcgen(
    inputSignals, 
    nResamplePoints=200, 
    CorridorRes=100)

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('Displacement (mm)')
plt.ylabel('Force (N)')
plt.title('Capsular Ligament')
fig.savefig('outputs/Corridors - Capsular Ligament.png')

# Interspinous Ligament
inputPath = 'Mattucci_CervicalLigaments/Mattucci_InterspinousLigament_QuasiStatic.csv'
rawCsv = np.genfromtxt(inputPath, delimiter=',', encoding=None)
nSignals = rawCsv.shape[1]/2
inputSignals = []
for iSig in range(int(nSignals)):
    temp = rawCsv[:,(2*iSig):(2*iSig+2)]
    index = np.logical_not(np.isnan(temp[:,0]))
    inputSignals.append(temp[index,:])

charAvg, innCorr, outCorr, processed, debug = arcgen(
    inputSignals, 
    nResamplePoints=200, 
    CorridorRes=100)

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('Displacement (mm)')
plt.ylabel('Force (N)')
plt.title('Interspinous Ligament')
fig.savefig('outputs/Corridors - Interspinous Ligament.png')

# Ligamentum Flavum
inputPath = 'Mattucci_CervicalLigaments/Mattucci_LigamentumFlavum_QuasiStatic.csv'
rawCsv = np.genfromtxt(inputPath, delimiter=',', encoding=None)
nSignals = rawCsv.shape[1]/2
inputSignals = []
for iSig in range(int(nSignals)):
    temp = rawCsv[:,(2*iSig):(2*iSig+2)]
    index = np.logical_not(np.isnan(temp[:,0]))
    inputSignals.append(temp[index,:])

charAvg, innCorr, outCorr, processed, debug = arcgen(
    inputSignals, 
    nResamplePoints=200, 
    CorridorRes=100)

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('Displacement (mm)')
plt.ylabel('Force (N)')
plt.title('Ligamentum Flavum')
fig.savefig('outputs/Corridors - Ligamentum Flavum.png')