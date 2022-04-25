from arcgen import arcgen
import numpy as np
from pathlib import Path

# inputPath = 'ExampleCasesAndDatasets/NBDL_15gFrontalDeceleration/NBDL_15gFrontal_HeadZAccel.csv'
inputPath = 'NBDL_15gFrontalDeceleration/NBDL_15gFrontal_HeadZAccel.csv'

avg, inn, out, processed, debug = arcgen(
    inputPath, 
    nResamplePoints=100, 
    CorridorRes=100, 
    nWarpCtrlPts=1,
    resultsToFile=False,
    Diagnostics='off')

rawCsv = np.genfromtxt(inputPath, delimiter=',', encoding=None)
nSignals = rawCsv.shape[1]/2
inputSignals = []
for iSig in range(int(nSignals)):
    temp = rawCsv[:,(2*iSig):(2*iSig+2)]
    index = np.logical_not(np.isnan(temp[:,0]))
    inputSignals.append(temp[index,:])

avg, inn, out, processed, debug = arcgen(
    inputSignals, 
    nResamplePoints=100, 
    CorridorRes=100, 
    nWarpCtrlPts=1,
    resultsToFile=False,
    Diagnostics='off')