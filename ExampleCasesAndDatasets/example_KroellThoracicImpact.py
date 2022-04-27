"""
Example: Kroell Thoracic Impact

This test cases demonstrates the effacacy of ARCGen on hysteretic or 
load-unload signals. In this case, the number of warping popints must be
chosen by examining the (x, s^hat) and (y, s^hat) for critical features
and inflection points. 

Please refer to test case readme document for more details. 

Dataset Citation:
    Kroell, C. K., Schneider, D. C., & Nahum, A. M. (1971). "Impact 
       Tolerance and Response of the Human Thorax." SAE Technical Papers.
       https://doi.org/10.4271/710851</div>

    Lobdell, T. E., Kroell, C. K., Schneider, D. C., Hering, W. E., & 
       Hahum, A. M. (1972). "Impact Response of the Human Thorax." In 
       W. F. King & H. J. Mertz (Eds.), "Human Impact Response: 
      Measurement and Simulation" (pp. 201–245). Springer Science + 
      Business Media.

Copyright (c) 2022 Devon C. Hartlen

"""

from arcgen import arcgen
import numpy as np
import matplotlib.pyplot as plt
import os

if not os.path.exists('outputs'):
    os.makedirs('outputs')


inputPath = 'Kroell_ThoracicImpactResponse/KroellThoracicImpactResponse.csv'
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
    Diagnostics='on',
    nWarpCtrlPts=2,
    WarpingPenalty=1e-2)

fig = plt.figure(figsize=(6,4), dpi=1200)
for signal in processed:
    pSig, = plt.plot(signal['resampled'][:,1], signal['resampled'][:,2], 
    label= 'Input Signals', c = 'lightgrey', lw=0.5)
pAvg, = plt.plot(charAvg[:,0], charAvg[:,1], label='Average', c ='black', lw = 1.5)
pInn, = plt.plot(innCorr[:,0], innCorr[:,1], label='Inner Corridor', c ='gold', lw = 1.5)
pOut, = plt.plot(outCorr[:,0], outCorr[:,1], label='Outer Corridor', c ='goldenrod', lw = 1.5)
fig.legend(handles=[pSig, pAvg, pInn, pOut])
plt.xlabel('Displacement (in)')
plt.ylabel('Force (lb)')
plt.title('Kroell Thoracic Impact Corridors')
fig.savefig('outputs/Corridors - Kroell Thoracic Impact.png')