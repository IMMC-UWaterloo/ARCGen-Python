from arcgen import arcgen
from pathlib import Path

inputPath = 'NBDL_15gFrontalDeceleration/NBDL_15gFrontal_HeadZAccel.csv'

avg, inn, out, processed, debug = arcgen(
    inputPath, 
    nResamplePoints=100, 
    CorridorRes=100, 
    nWarpCtrlPts=1,
    resultsToFile=False,
    Diagnostics='off')

