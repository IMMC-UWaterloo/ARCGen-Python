# %% ARCGen - Python | Arc-length Response Corridor Generator
# %
# % Created By:     D.C. Hartlen, M.ASc, EIT
# % Date:           27-Jun-2021
# % Updated By:     D.C. Hartlen, M.ASc, EIT
# % Date:           24-Feb-2022
# % Version:        Python 3.x
# %
# % ARCGen, short for Arc-length Response Corridor Generation, provides
# % automated calculation of a characteristic average and response corridors
# % on input signals regardless of if said signals are non-monotonic or
# % hystertic. This is accomplished by re-parameterizing input signals based
# % on arc-length. Corridors are extracted using a marching squares
# % algorithm.
# %
# % If you use ARCGen in your research, please use the following citation:
# %     Hartlen, D.C. & Cronin, D.S. (2022). "Arc-length Re-parametrization
# %        and Signal Registration to Determine a Characteristic Average and
# %        Statistical Response Corridors of Biomechanical Data." Frontiers
# %        in Bioengineering and Biotechnology.
# %
# % ARCGen is released under a GNU GPL v3 license. No warranty or support is
# % provided. The authors hold no responsibility for the validity, accuracy,
# % or applicability of any results obtained from this code.
# %
# % This function has one mandatory input, four outputs, and many optional
# % inputs. Optional inputs are defined using name-value pair arguments.
# %
# % Usage notes:
# % It is common to see errors when running this function if the number of
# % resampling points or corridor resolution is too sparse or signals
# % exhibit significant variablity not accounted for through signal
# % registration. This tends to manifest in either truncated corridors or the
# % code termininating in an error. Often increasing resampling points or
# % corridor resolution. Turning 'Diagnostics' to 'detailed' can help
# % identify these issues.
# %
# % Computed corridors will often not extend all the way to the shared origin
# % of input signals. This is because small low st. dev. at this shared point
# % is too low to be captured during corridor extraction with the marching
# % squares algorithm. There are two solutions to this problem. First, one
# % could force a minimum corridors size using the 'MinCorridorWidth' option.
# % Second, one could manually extend corridors in post-processing.
# %
# % MANDATORY INPUTS:
# % -----------------
# % A csv file including the signals to be placed in the working directory. The signals (data) needs to be placed in the csv file with or with a header
#   such that each signal has a column for x-axis data and a column for y-axis data.
# %
# % OPTIONAL INPUTS:
# % ----------------
# % nResamplePoints: integer defining the number of points used to
# %       re-parameterize input signals. Default: 100.
# % CorridorRes: integeer defining the number of grid points used for the
# %       marching squares algorithm. The sampling grid for the marching
# %       squares algorithm extends 120% of extreme corridors. This parameter
# %       defines the number of points along each side of the grid.
# %       Default: 100. It is common to increase this significantly.
# % NormalizeSignals: character arry used to turn on signal normalization.
# %       Options: 'on' (default), 'off'
# % EllipseKFact: float used to scale the major and minor axis of the
# %       ellipses used for corridor generation. This value corrisponds to
# %       the square root of the chi-squared CDF. Default: 1.0 (creates
# %       corridors one standard deviation along the x and y axes)
# % Diagnostics: character array used to activate diagnostic plots. Useful
# %       for debugging errors. Options: 'off' (default), 'on', 'detailed'.
# % MinCorridorWidth: Factor used to enforce a minimum corridor width. Any
# %       st.dev. less than 'MinCorridorFactor'*max(st.dev.) is replaced with
# %       'MinCorridorFactor'*max(st.dev.). x & y axes are handled
# %       separately. A value of 0 (default) disables forcing minimum width.
# % nWarpCtrlPts: integer that sets the number of interior control points
# %       used for signal registration. A value of 0 (default) disables
# %       signal registration
# % WarpingPenalty: float specifying the penalty factor used during the
# %       signal registration process. A value of 10^-2 (default) to 10^3 is
# %       recommended, but the exact value will need to be tuned to a
# %       specific problem.
# %
# % Mandatory OUTPUTS:
# % ------------------
# % charAvg: an [nResamplePoints,2] array containing the computed
# %       characteristic average.
# % innerCorr: an [nResamplePoints,2] array containing points defining the
# %       inner corridor
# % outerCorr: an [nResamplePoints,2] array containing points defining the
# %       outer corridor
# % Plot: a plot of the input signal, characteristic average and signal corridors.


## Modules
from statistics import stdev
import numpy as np
from numpy import inf, inner
import arcgen.polygonFunctions as poly
from arcgen.uniquetol import uniquetol
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy import interpolate
from scipy import optimize
import os
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

def arcgen(inputDataPath,
           nResamplePoints = 250,
           Diagnostics = 'off',
           NormalizeSignals = 'on',
           EllipseKFact = 1,
           CorridorRes = 250,
           MinCorridorWidth = 0,
           nWarpCtrlPts = 4,
           WarpingPenalty = 1e-2,
           resultsToFile = False):

  # Creating a directory for results
  if (Diagnostics == 'on') or (Diagnostics == 'detailed') or (resultsToFile):
    if not os.path.exists('outputs'):
      os.makedirs('outputs')

  # Declarations, as dictionaries
  inputSignals ={}
  maxAlen = {}
  meanDev = {}
  medianDevs = {}
  normalizedSignal = {}
  xMax = {}
  xMin = {}
  yMax = {}
  yMin = {}
  xNormMax = {}
  yNormMax = {}
  warpControlPoints = {}

  if Path(inputDataPath).exists():
    # Load from csv. empty values in shorter signals are given as 'nan'
    dataframe = np.genfromtxt(inputDataPath, delimiter=',', encoding=None)
    numberRows, numberCols = dataframe.shape
    # Error check
    if numberCols % 2 == 0:
      numberSignals = int(numberCols/2)
    else:
      raise ValueError("The number of columns in the csv file is not even")    
  else:
    raise ValueError("inputData format is not valid. Only valid paths or tuples") 

  # Store input signals as list of arrays
  for i in range(len(np.hsplit(dataframe,numberSignals))):
    temp = dataframe[:, (2*i):(2*i+2)]
    # Remove 'nan' entries
    indexNan = np.logical_not(np.isnan(temp[:,0]))
    # add to dictionary
    inputSignals['Signal '+str(i+1)] = temp[indexNan,:]

  #%% Compute arclength based on input signal datapoints
  #% Do not perform normalization
  if NormalizeSignals == 'off':
    for iSignal in inputSignals.keys():
      temp = inputSignals[iSignal].copy() # Temporary for convenience

      # Compute arc - length between each data point
      segments = np.sqrt((temp[0:-1,0] - temp[1:,0])**2 + (temp[0:-1, 1] - temp[1:, 1])**2)
      segments = np.append([0], segments)

      #% Append cumulative arc length to data array
      alen = np.cumsum(segments)

      #% Compute normalized arc - length
      maxAlen[iSignal] = np.max(alen)
      inputSignals[iSignal] = np.column_stack((inputSignals[iSignal], alen))
      inputSignals[iSignal] = np.column_stack((inputSignals[iSignal], alen/maxAlen[iSignal]))

      #% Determine max[x, y] data
      xMax[iSignal], yMax[iSignal] = temp.max(axis=0)

      #% Remove spurious duplicates
      _,idx=np.unique(inputSignals[iSignal], axis=0,return_index=True)
      inputSignals[iSignal][np.sort(idx)]

  #% Perform magnitude normalization based on bounding box
  if NormalizeSignals =='on':
    #% Determine bounding box of individual signals
    for iSignal in inputSignals.keys():
      temp = inputSignals[iSignal]  # Temporary for convenience

  #% Determine min[x, y] and max[x, y]data
      xMin[iSignal], yMin[iSignal] = temp.min(axis=0)
      xMax[iSignal], yMax[iSignal] = temp.max(axis=0)

    xBound = [sum(xMin.values())/len(xMin.values()), sum(xMax.values())/len(xMax.values())]
    yBound = [sum(yMin.values())/len(yMin.values()), sum(yMax.values())/len(yMax.values())]

    # Normalize the axis of each signal, then do arc-length calcs
    for iSignal in inputSignals.keys():
      # This needs to be a .copy(), otherwise scaling effects inputSignals 
      temp = inputSignals[iSignal].copy() # Temporary for convenience

      #% Normalize from bounding box to [-1,1]
      temp[:,0] = temp[:,0] / (xBound[1]-xBound[0])
      temp[:,1] = temp[:,1] / (yBound[1]-yBound[0])

      # % Compute arc - length between each data point
      segments = np.sqrt((temp[0:-1, 0] - temp[1:, 0]) ** 2 + (temp[0:-1, 1] - temp[1:, 1])** 2)
      segments = np.append([0], segments)

      # % Append cumulative arc length to data array
      alen = np.cumsum(segments)

      # % Compute normalized arc - length
      maxAlen[iSignal] = np.max(alen)
      inputSignals[iSignal] = np.column_stack((inputSignals[iSignal], alen))
      inputSignals[iSignal] = np.column_stack((inputSignals[iSignal], alen / maxAlen[iSignal]))

      # % Determine max[x, y] data
      xNormMax[iSignal], yNormMax[iSignal] = temp.max(axis=0)

      # % Remove spurious duplicates
      temp = inputSignals[iSignal].copy()
      _, uniqueIndex = np.unique(temp[:, 3], return_index=True)
      inputSignals[iSignal] = temp[uniqueIndex, :]

  #% Compute mean and median arc - length deviation
  meanAlen = sum(maxAlen.values())/len(maxAlen.values())
  medianAlen = np.median(np.fromiter(maxAlen.values(), dtype=float))

  for iSignal in inputSignals.keys():
    meanDev[iSignal] = maxAlen[iSignal] - meanAlen
    medianDevs[iSignal] = maxAlen[iSignal] - medianAlen

    normAlen = np.linspace(0, inputSignals[iSignal][-1, 3], num = nResamplePoints)

    resampX = interpolate.interp1d(inputSignals[iSignal][:, 3], inputSignals[iSignal][:, 0])(normAlen)
    resampY = interpolate.interp1d(inputSignals[iSignal][:, 3], inputSignals[iSignal][:, 1])(normAlen)

    #% Resulting array is normalized arc - length, resampled x, resam.y
    normalizedSignal[iSignal] = np.column_stack((normAlen, resampX, resampY))

  #% % For each resampled point, determine average and standard deviation across signals
  #% Initialize arrays
  charAvg = np.zeros((nResamplePoints, 2))
  stdevData = np.zeros((nResamplePoints, 2))

  for iPoints in range(nResamplePoints):
    temp = np.empty([0, 2])

    for iSignal in inputSignals.keys():
      temp = np.vstack((temp, normalizedSignal[iSignal][iPoints,1:3]))
    charAvg[iPoints,:]= temp.mean(axis=0)
    stdevData[iPoints,:] = temp.std(axis=0)

  #%% Align normalized arc-length signals based on minimized correlation.
  #% Enabled by option 'nWarpCtrlPts'. If 0, skip alignment.
  if nWarpCtrlPts > 0:
    #% Assemble signal matrices prior to correlation
    signalX = np.zeros((nResamplePoints, len(inputSignals)))
    signalY = np.zeros((nResamplePoints, len(inputSignals)))
    for i in range(len(inputSignals)):
      signalX[:,i] = normalizedSignal[list(inputSignals)[i]][:, 1]
      signalY[:,i] = normalizedSignal[list(inputSignals)[i]][:, 2]
    meanCorrScore, corrArray = evalCorrScore(signalX,signalY)
    #% Assign pre-optimized correlation scores to debug structure
    preWarpCorrArray = corrArray;
    preWarpMeanCorrScore = meanCorrScore;

    #% Optimize warp points for arbitrary n warping points. Build bounds,
    #% constraints, and x0s
    nWarp = nWarpCtrlPts
    nSignal = numberSignals

    if nWarp == 1:   #% nWarp == 1 is a special case as inequalites aren't needed
      x0 = 0.5 * np.ones((nSignal*2))
      lb = 0.15 * np.ones((nSignal*2))
      ub = 0.85 * np.ones((nSignal*2))
      A = []
      b = []
      linConstraints = []
    elif nWarp >= 15:
      raise ValueError('Specifying more than 10 interior warping points is not supported')
    else:
      x0 = np.zeros((nWarp * (nSignal * 2)))

      for i in range(nWarp):
        x0[np.array([((i) * nSignal)+ np.arange(0, nSignal)+(i * (nSignal))])] = ((i+1)/(nWarp+1)) * np.ones((nSignal))
        x0[np.array([((i) * nSignal) + np.arange(0, nSignal) + ((i+1) * nSignal)])] = ((i+1)/(nWarp+1)) * np.ones((nSignal))

      lb = 0.05 * np.ones((nWarp*(nSignal*2)))
      ub = 0.95 * np.ones((nWarp*(nSignal*2)))
      A = np.zeros([((nWarp-1) * (nSignal * 2)), (nWarp * (nSignal * 2))]) # may not be needed but added for future use
      b = -0.05 * np.ones(((nWarp-1) * (nSignal * 2)))  # % Force some separation between warped points
      
      for iSignal in range(nSignal * 2):
        for iWarp in range(nWarp-1):
          A[iSignal + (iWarp) * (nSignal * 2), iSignal + ((iWarp) * (nSignal * 2))] = 1
          A[iSignal + (iWarp) * (nSignal * 2), iSignal + ((iWarp+1) * (nSignal * 2))] = -1
      
      linConstraints = optimize.LinearConstraint(A, -np.inf *np.ones(((nWarp-1) * (nSignal * 2))), b)
    
    #% Execute optimization and compute warped signals
    bounds = list(zip(lb, ub))
    optWarpArrayX = optimize.minimize(warpingObjective, x0, args=(nWarp,inputSignals, WarpingPenalty,nResamplePoints), method = 'SLSQP', bounds=bounds, options={ 'disp': False}, constraints=linConstraints)
    optWarpArray =optWarpArrayX.x
    optWarpArray = optWarpArray.reshape((nSignal*2,nWarp), order='F')

    warpedSignals, signalX, signalY = warpArcLength(optWarpArray,inputSignals,nResamplePoints)


    #% Compute correlation score
    meanCorrScore, corrArray = evalCorrScore(signalX, signalY)
    #% Assign warped correlation scores to debug structure
    warpedCorrArray = corrArray
    warpedMeanCorrScore = meanCorrScore

        
    #% Replace 'normalizedSignal' in 'responseSignal' and compute average and
    #% standard deviation.
    for iSignal, keys in enumerate(inputSignals):
      normalizedSignal[keys] = warpedSignals[keys]
      tempRow1 = np.concatenate([[0],optWarpArray[iSignal+nSignal,:],[1]])
      tempRow2 = np.concatenate([[0],optWarpArray[iSignal,:],[1]])
      warpControlPoints[keys] = np.concatenate([[tempRow1], [tempRow2]])
      
    for iPoints in range(nResamplePoints):
      temp = np.empty((nSignal, 2)) #% probably cleaner way to do this.
      for iSignal, keys in enumerate(inputSignals):
        #% collect specific point from each signal
        temp[iSignal,:] = normalizedSignal[keys][iPoints,1:3]
      charAvg[iPoints,:] = np.mean(temp, axis = 0)
      stdevData[iPoints,:] = np.std(temp, axis = 0,  ddof=1)

  #%% Clamp minimum corridor width. Disabled if 'MinCorridorWidth' == 0
  #% Include influence of corridor scaling factor 'EllipseKFact'
  if MinCorridorWidth > 0:
    #% Replace any stDevData below maximum st.dev. * 'MinCorridorWidth'
    stdevData[:,0] = np.where(stdevData[:,0] < (MinCorridorWidth * np.max(stdevData[:,0]) * EllipseKFact), (MinCorridorWidth * np.max(stdevData[:,0]) * EllipseKFact), stdevData[:,0])
    stdevData[:,1] = np.where(stdevData[:,1] < (MinCorridorWidth * np.max(stdevData[:,1]) * EllipseKFact), (MinCorridorWidth * np.max(stdevData[:,1]) * EllipseKFact), stdevData[:,1])

  #%% Diagnostic: Plot normalized signals and St. Devs. 
  if Diagnostics == 'on' or Diagnostics == 'detailed':
      
    #% Plot normalized x,y 
    fig, ax = plt.subplots(2, 2, figsize = (8,6), dpi = 300)
    fig.suptitle('Arc-length Discretized Normalized Signals')
    ax[0][0].title.set_text('Arc-length Discretized Normalized Signals')
    ax[0][1].title.set_text('Warping Functions')
    ax[1][0].title.set_text('Average and St.Dev. of X-Data')
    ax[1][1].title.set_text('Average and St.Dev. of Y-Data')


    for iSignal, keys in enumerate(inputSignals):
      ax[0][0].plot(normalizedSignal[keys][:,1], normalizedSignal[keys][:,2], label= keys)
      
      if nWarpCtrlPts > 0:
        interpX = np.concatenate([[0],optWarpArray[iSignal+nSignal,:],[1]], axis=None, dtype='float')
        interpY = np.concatenate([[0],optWarpArray[iSignal,:],[1]], axis=None, dtype='float')
        ax[0][1].plot(interpX, interpY, lw=0.0, marker='x', markersize=6)
        interpResults = interpolate.pchip(interpX, interpY, axis=0)(inputSignals[keys][:,3])
        ax[0][1].plot(inputSignals[keys][:,3],interpResults, label = keys )
      else: 
          ax[0][1].title.set_text('No Warping Performed')

      ax[1][0].plot(normalizedSignal[list(normalizedSignal.keys())[0]][:,0], normalizedSignal[keys][:,1], label= keys)
      ax[1][0].errorbar(normalizedSignal[list(normalizedSignal.keys())[0]][:,0], charAvg[:,0], yerr=stdevData[:,0])
      ax[1][1].plot(normalizedSignal[list(normalizedSignal.keys())[0]][:,0], normalizedSignal[keys][:,2], label= keys)
      ax[1][1].errorbar(normalizedSignal[list(normalizedSignal.keys())[0]][:,0], charAvg[:,1], yerr=stdevData[:,1])


      
    #ax[0][1].plot(interpX, interpY )
    ax[0][1].plot([0,1],[0,1])
    ax[0][0].set(xlabel='x-data', ylabel='y-data')
    ax[0][1].set(xlabel='Unwarped Normalized Arc-length', ylabel='Warped Normalized Arc-length')
    ax[1][0].set(xlabel='Normalized Arc-length', ylabel='x-data')
    ax[1][1].set(xlabel='Normalized Arc-length', ylabel='y-data')



    fig.tight_layout()
    fig.subplots_adjust(top=0.91)
    ax[0][0].legend(loc='lower right')
    fig.savefig('outputs/Normalization and Warping.png')

  #%% Begin marching squares algorithm
  #% Create grids based on upper and lower of characteristic average plus 120%
  #% of maximum standard deviation
  scaleFact = 1.25 * EllipseKFact
  xx,yy = np.meshgrid(np.linspace((np.min(charAvg[:,0]) - scaleFact*np.max(stdevData[:,0])), (np.max(charAvg[:,0]) + scaleFact * np.max(stdevData[:,0])), num = CorridorRes), np.linspace((np.min(charAvg[:,1]) - scaleFact*np.max(stdevData[:,1])), (np.max(charAvg[:,1]) + scaleFact*np.max(stdevData[:,1])), num = CorridorRes), indexing='xy')
  zz = np.zeros(np.shape(xx))   #% initalize grid of ellipse values
  # #% For each grid point, find the max of each standard deviation ellipse
  kFact = EllipseKFact #% faster if no struct call in inner loop. 
  nRes = CorridorRes   #% again, for speed
  for iPt in range(nRes):
    for jPt in range(nRes):
      zz[iPt,jPt] = np.max(
          (((xx[iPt,jPt] - charAvg[:,0])**2 / (stdevData[:,0]*kFact)**2 + (yy[iPt,jPt] - charAvg[:,1])**2 / (stdevData[:,1]*kFact)**2)**-1))
      
  # % The following segments is the marching squares algorith. The goal of this
  # % algorithm is to find the zz=1 isoline, as this represents the outer
  # % boundary of all elllipses. 
  # %
  # % Described in brief, this algorithm goes through each point, looking at
  # % its and its neighbours values. There are only 16 configurations of these
  # % squares or cells. Based on the configuration, add the appropriate line
  # % segments. This method uses linear interpolation to increase accuracy. 
  # % Initalize line segments for speed. This line may cause issues, as it
  # % assumes maximum size. Bump up 10 if it does. 
  lineSegments = np.zeros((10 * max(nResamplePoints,CorridorRes), 4))

  iSeg = -1

  for iPt in range(CorridorRes-1):  #% Rows (y-axis)
    for jPt in range(CorridorRes-1):   #% Columns (x-axis)
      # % Cell value definition
      # %  1 -- 2 
      # %  |    |
      # %  |    |
      # %  8 -- 4
      # %
      # % REMEMBER!!!! 
      # % array(i,j) = array(rows, columns,) = array(y,x)
      
      # % By carefully defining cell values and definitions, we can use
      # % binary to simplify logic though a integer based switch case

      cellValue = int((1 * (zz[iPt,jPt]>1)) + (2 * (zz[iPt+1,jPt]>1)) + (4 * (zz[iPt+1,jPt+1]>1)) + (8 * (zz[iPt,jPt+1]>1)) + 1)
      
      if cellValue == 1:
        #% No Vertices
        pass

      elif  cellValue == 2:
        #% South-West
        iSeg = iSeg+1
        lineSegments[iSeg,:] = [interpIso(xx[iPt,jPt],zz[iPt,jPt],xx[iPt,jPt+1],zz[iPt,jPt+1]),yy[iPt,jPt],
        xx[iPt,jPt], interpIso(yy[iPt,jPt],zz[iPt,jPt],yy[iPt+1,jPt],zz[iPt+1,jPt])]

      elif  cellValue == 3:
        #% West-North
        iSeg = iSeg+1
        lineSegments[iSeg,:] = [xx[iPt+1,jPt],interpIso(yy[iPt,jPt],zz[iPt,jPt], yy[iPt+1,jPt],zz[iPt+1,jPt]),
            interpIso(xx[iPt+1,jPt],zz[iPt+1,jPt], xx[iPt+1,jPt+1],zz[iPt+1,jPt+1]),yy[iPt+1,jPt]]
      
      elif  cellValue == 4:
        #% North-South
        iSeg = iSeg+1;
        lineSegments[iSeg,:] = [interpIso(xx[iPt,jPt],zz[iPt,jPt],xx[iPt,jPt+1],zz[iPt,jPt+1]),yy[iPt,jPt],
            interpIso(xx[iPt+1,jPt],zz[iPt+1,jPt],xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),yy[iPt+1,jPt]]
            
      elif  cellValue == 5:
        #% North-East
        iSeg = iSeg+1;
        lineSegments[iSeg,:] = [interpIso(xx[iPt+1,jPt],zz[iPt+1,jPt],xx[iPt+1,jPt+1],zz[iPt+1,jPt+1]),yy[iPt+1,jPt+1],
            xx[iPt+1,jPt+1],interpIso(yy[iPt+1,jPt+1],zz[iPt+1,jPt+1], yy[iPt,jPt+1], zz[iPt,jPt+1])]
            
      elif  cellValue == 6:  #% Ambiguous 
        centerVal = np.mean(list([zz[iPt,jPt], zz[iPt+1,jPt], zz[iPt+1,jPt+1], zz[iPt, jPt+1]]))
        if centerVal >= 1:
          #% West-North
          iSeg = iSeg+1
          lineSegments[iSeg,:] = [xx[iPt+1,jPt],interpIso(yy[iPt,jPt],zz[iPt,jPt], yy[iPt+1,jPt],zz[iPt+1,jPt]),
              interpIso(xx[iPt+1,jPt],zz[iPt+1,jPt], xx[iPt+1,jPt+1],zz[iPt+1,jPt+1]),yy[iPt+1,jPt]]
                
          #% South-East
          iSeg = iSeg+1
          lineSegments[iSeg,:] = [interpIso(xx[iPt,jPt+1],zz[iPt,jPt+1],xx[iPt,jPt],zz[iPt,jPt]),yy[iPt,jPt+1],
              xx[iPt,jPt+1],interpIso(yy[iPt,jPt+1],zz[iPt,jPt+1],yy[iPt+1,jPt+1],zz[iPt+1,jPt+1])]
        else:
          #% South-West
          iSeg = iSeg+1
          lineSegments[iSeg,:] = [interpIso(xx[iPt,jPt],zz[iPt,jPt],xx[iPt,jPt+1],zz[iPt,jPt+1]),yy[iPt,jPt],
              xx[iPt,jPt], interpIso(yy[iPt,jPt],zz[iPt,jPt],yy[iPt+1,jPt],zz[iPt+1,jPt])]
                
          #% North-East
          iSeg = iSeg+1
          lineSegments[iSeg,:] = [interpIso(xx[iPt+1,jPt],zz[iPt+1,jPt],xx[iPt+1,jPt+1],zz[iPt+1,jPt+1]),yy[iPt+1,jPt+1],
              xx[iPt+1,jPt+1],interpIso(yy[iPt+1,jPt+1],zz[iPt+1,jPt+1], yy[iPt,jPt+1], zz[iPt,jPt+1])]

      elif  cellValue == 7:
        #% West-East
        iSeg = iSeg+1
        lineSegments[iSeg,:] = [xx[iPt,jPt],interpIso(yy[iPt,jPt],zz[iPt,jPt],yy[iPt+1,jPt],zz[iPt+1,jPt]),
            xx[iPt,jPt+1],interpIso(yy[iPt,jPt+1],zz[iPt,jPt+1],yy[iPt+1,jPt+1],zz[iPt+1,jPt+1])]
      
      elif  cellValue == 8:
        #% South - East
        iSeg = iSeg+1
        lineSegments[iSeg,:] = [interpIso(xx[iPt,jPt+1],zz[iPt,jPt+1],xx[iPt,jPt],zz[iPt,jPt]),yy[iPt,jPt+1],
            xx[iPt,jPt+1],interpIso(yy[iPt,jPt+1],zz[iPt,jPt+1],yy[iPt+1,jPt+1],zz[iPt+1,jPt+1])]
                                      
      elif  cellValue == 9:
        #% South - East
        iSeg = iSeg+1
        lineSegments[iSeg,:] = [interpIso(xx[iPt,jPt+1],zz[iPt,jPt+1],xx[iPt,jPt],zz[iPt,jPt]),yy[iPt,jPt+1],
            xx[iPt,jPt+1],interpIso(yy[iPt,jPt+1],zz[iPt,jPt+1],yy[iPt+1,jPt+1],zz[iPt+1,jPt+1])]

      elif  cellValue == 10:
        #% West-East
        iSeg = iSeg+1
        lineSegments[iSeg,:] = [xx[iPt,jPt],interpIso(yy[iPt,jPt],zz[iPt,jPt],yy[iPt+1,jPt],
            zz[iPt+1,jPt]), xx[iPt,jPt+1],interpIso(yy[iPt,jPt+1],zz[iPt,jPt+1],yy[iPt+1,jPt+1],zz[iPt+1,jPt+1])]

      elif  cellValue == 11: #% Ambiguous
        centerVal = np.mean(list([zz[iPt,jPt], zz[iPt+1,jPt], zz[iPt+1,jPt+1], zz[iPt, jPt+1]]))
        if centerVal >= 1:
          #% South-West
          iSeg = iSeg+1
          lineSegments[iSeg,:] = [interpIso(xx[iPt,jPt],zz[iPt,jPt],xx[iPt,jPt+1],zz[iPt,jPt+1]),yy[iPt,jPt],
              xx[iPt,jPt], interpIso(yy[iPt,jPt],zz[iPt,jPt],yy[iPt+1,jPt],zz[iPt+1,jPt])]

          #% North-East
          iSeg = iSeg+1
          lineSegments[iSeg,:] = [interpIso(xx[iPt+1,jPt],zz[iPt+1,jPt],xx[iPt+1,jPt+1],zz[iPt+1,jPt+1]),yy[iPt+1,jPt+1],
              xx[iPt+1,jPt+1],interpIso(yy[iPt+1,jPt+1],zz[iPt+1,jPt+1], yy[iPt,jPt+1], zz[iPt,jPt+1])]
        else:
          #% West-North
          iSeg = iSeg+1
          lineSegments[iSeg,:] = [xx[iPt+1,jPt],interpIso(yy[iPt,jPt],zz[iPt,jPt], yy[iPt+1,jPt],zz[iPt+1,jPt]),
                interpIso(xx[iPt+1,jPt],zz[iPt+1,jPt], xx[iPt+1,jPt+1],zz[iPt+1,jPt+1]),yy[iPt+1,jPt]]
          #% South-East
          iSeg = iSeg+1
          lineSegments[iSeg,:] = [interpIso(xx[iPt,jPt+1],zz[iPt,jPt+1],xx[iPt,jPt],zz[iPt,jPt]),yy[iPt,jPt+1],
              xx[iPt,jPt+1],interpIso(yy[iPt,jPt+1],zz[iPt,jPt+1],yy[iPt+1,jPt+1],zz[iPt+1,jPt+1])]

      elif  cellValue == 12:
        #% North-East
        iSeg = iSeg+1
        lineSegments[iSeg,:] = [interpIso(xx[iPt+1,jPt],zz[iPt+1,jPt],xx[iPt+1,jPt+1],zz[iPt+1,jPt+1]),yy[iPt+1,jPt+1],
            xx[iPt+1,jPt+1],interpIso(yy[iPt+1,jPt+1],zz[iPt+1,jPt+1], yy[iPt,jPt+1], zz[iPt,jPt+1])]
      
      elif  cellValue == 13:
        #% North-South
        iSeg = iSeg+1
        lineSegments[iSeg,:] = [interpIso(xx[iPt,jPt],zz[iPt,jPt],xx[iPt,jPt+1],zz[iPt,jPt+1]),yy[iPt,jPt],
            interpIso(xx[iPt+1,jPt],zz[iPt+1,jPt],xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),yy[iPt+1,jPt]]
      
      elif  cellValue == 14:
        #% West-North
        iSeg = iSeg+1
        lineSegments[iSeg,:] = [xx[iPt+1,jPt],interpIso(yy[iPt,jPt],zz[iPt,jPt], yy[iPt+1,jPt],zz[iPt+1,jPt]),
            interpIso(xx[iPt+1,jPt],zz[iPt+1,jPt], xx[iPt+1,jPt+1],zz[iPt+1,jPt+1]),yy[iPt+1,jPt]]

      elif  cellValue == 15:
        #% South-West
        iSeg = iSeg+1
        lineSegments[iSeg,:] = [interpIso(xx[iPt,jPt],zz[iPt,jPt],xx[iPt,jPt+1],zz[iPt,jPt+1]),yy[iPt,jPt],
            xx[iPt,jPt], interpIso(yy[iPt,jPt],zz[iPt,jPt],yy[iPt+1,jPt],zz[iPt+1,jPt])]

      elif  cellValue == 16:
        pass
        #% No vertices

  lineSegments = lineSegments[0:iSeg+1,:]

  #% Extract list of unique vertices from line segmens
  vertices = np.vstack([lineSegments[:,:2],lineSegments[:,2:4]])
  vertices = uniquetol(vertices, 1e-12)
  index = np.zeros(len(vertices), dtype=bool)

  #% Create a vertex connectivity table. The 1e-12 value is here because
  #% floats don't round well and == will not work. 
  vertConn = np.zeros((len(lineSegments),2)).astype(int)

  for i in range(len(vertConn)):
    index = np.all(np.abs(lineSegments[:,:2] - vertices[i,:]) < 1e-12, axis = 1)
    vertConn[index, 0] = i
    index = np.all(np.abs(lineSegments[:,2:4] - vertices[i,:]) < 1e-12, axis = 1)
    vertConn[index, 1] = i

  #%% Start line segments sorting and envelope extraction
  nEnvelopes = 0
  allEnvelopes = np.zeros((len(vertConn),2))

  for i in range(len(vertConn)-1):
    j = i + 1 #% helper index

    #% save vertex to find 
    vertToFind = vertConn[i,1].astype(int)

    #% Find connecting node
    vertConnToFind_ = np.any((vertConn[j:,:].astype(int) == vertToFind.astype(int)), axis = 1)
    foundShiftedInd = find(vertConnToFind_.astype(int), 'first')

    #% If we have found an index 
    if foundShiftedInd == 0 and vertConn[j,0] == vertToFind:
      continue

    elif foundShiftedInd == 0 and vertConn[j,1] == vertToFind:
      vertConn[j, [0, 1]] = vertConn[j, [1, 0]]
      continue
    
    elif foundShiftedInd > 0:
      foundInd = foundShiftedInd + j

      #% swap found vert conn row with j row
      vertConn[[j,foundInd],:] = vertConn[[foundInd,j],:]
      if vertConn[j,0] == vertToFind:
        continue
      else:
        vertConn[j,[0,1]] = vertConn[j, [1,0]]
      

    #% If we did not find an index, we either may have an open envelope or
    #% envelope may be convex and loops back on itself. 
    elif foundShiftedInd.size == 0:
      #% Check to see if we can find the first vertex in envelope
      #% appearing again (check for closure)
      vertConn_ = allEnvelopes[nEnvelopes,0].astype(int)
      vertToFind = vertConn[vertConn_,1] ####### Confirm column
      foundShiftedInd = find(np.any(vertConn[j:,:] == vertToFind, axis = 1), 'first')
      
      #% If we do not find an index, it means this envelope is complete and manifold
      if foundShiftedInd.size == 0:
      #% Assign indices to finish current envelope, initialize next
        allEnvelopes[nEnvelopes,1] = i
        nEnvelopes = nEnvelopes + 1
        allEnvelopes[nEnvelopes, 0] = j
      else:
        raise ValueError('Literal Edge Case')

  allEnvelopes[nEnvelopes,1] = j

  #% Find largest envelope
  envInds = np.argmax((allEnvelopes[:,1]-allEnvelopes[:,0])).astype(int)
  #% Convert indices in evelopes to array of (x,y)
  envInds = allEnvelopes[envInds, :]
  envelope = vertices[vertConn[envInds[0].astype(int):envInds[1].astype(int)+1,0] ,:]

  #%% Divide the envelope into corridors. 
  #% To break the largest envelop into inner and outer corridors, we need to
  #% account for several edge cases. First, we test to see if there are any
  #% intercepts of the characteristic average and the largest envelope. 
  closedEnvelope = np.vstack((envelope, envelope[0,:]))

  _, indexIntercept = poly.polyxpoly(closedEnvelope,charAvg)

  indexIntercept = np.floor(indexIntercept).astype(int)

  #% If we find two intercepts, then we have no problem
  if indexIntercept.shape[0] >= 2:
    # Sort to order charAvg interecepts (2nd column), smaller charAvg intercepts
    # are closer to start. 
    indexSort = np.argsort(indexIntercept[:,1], axis=0)
    iIntStart = indexIntercept[indexSort[0],0]
    iIntEnd = indexIntercept[indexSort[-1],0]

  #% If we find only one intercept, we need to determine if the intercept is a
  #% the start or end of the envelope. Then we need to extend the opposite
  #% side of the characteristic average to intercept the envelope. 
  elif indexIntercept.shape[0] == 1:
    #% Compute extension 
    aLenInterval = 1/nResamplePoints
    indexLength = round(0.2*len(charAvg))
      
    aLenExtension = np.abs(aLenInterval/(charAvg[0,:]-charAvg[1,:])) * 1.1 * np.max(stdevData)
    aLenExtension[aLenExtension == inf] = 0
    aLenExtension = max(aLenExtension)
    
    
    #% If the single found point is inside the envelope, the found intercept
    #% is at the end. Therefore extend the start
    if poly.inpolygon(envelope, charAvg[indexIntercept[0,1],:] ): 
      iIntEnd = indexIntercept[-1,0]
      linestart_0 = interpLin(0, charAvg[0,0], aLenInterval, charAvg[1,0], -aLenExtension)
      linestart_1 = interpLin(0, charAvg[0,1], aLenInterval, charAvg[1,1], -aLenExtension)

      lineStart =  np.vstack((np.vstack(np.asarray([[linestart_0 , linestart_1]])), charAvg[0:indexLength,:]))

    #%Find intercepts to divide line using Poly
      _, iIntStart = poly.polyxpoly(closedEnvelope, lineStart)
      iIntStart = np.floor(iIntStart[0,0]).astype(int)

    #% If the single found point is outside the envelope, the found
    #% intercept is the start
    else:
      iIntStart = indexIntercept[0,0]
      # charAvg = np.asarray(charAvg)
      lineend_0 = interpLin(1-aLenInterval, charAvg[-2,0], 1, charAvg[-1,0], (1+aLenExtension))
      lineend_1 = interpLin(1-aLenInterval, charAvg[-2,1], 1, charAvg[-1,1], (1+aLenExtension))

      lineEnd =  np.vstack((charAvg[-1-indexLength:-1,:], np.vstack(np.asarray([[lineend_0 , lineend_1]]))))
      
      #%Find intercepts to divide line using Poly
      _, iIntEnd = poly.polyxpoly(closedEnvelope, lineEnd)
      iIntEnd = np.floor(iIntEnd[-1,0]).astype(int)
      
  #% If we find no intercepts, we need to extend both sides of characteristic
  #% average to intercept the envelop.
  else:
    aLenInterval = 1/nResamplePoints
    indexLength = round(0.2*len(charAvg))
    
    aLenExtension = np.abs(aLenInterval/(charAvg[0,:]-charAvg[1,:])) * 1.1 * np.max(stdevData)
    aLenExtension[aLenExtension == inf] = 0
    aLenExtension = max(aLenExtension)

    linestart_0 = interpLin(0, charAvg[0,0], aLenInterval, charAvg[1,0], -aLenExtension)
    linestart_1 = interpLin(0, charAvg[0,1], aLenInterval, charAvg[1,1], -aLenExtension)
    lineStart =  np.vstack((np.vstack(np.asarray([[linestart_0 , linestart_1]])), charAvg[0:indexLength,:]))
    
    lineend_0 = interpLin(1-aLenInterval, charAvg[-2,0], 1, charAvg[-1,0], (1+aLenExtension))
    lineend_1 = interpLin(1-aLenInterval, charAvg[-2,1], 1, charAvg[-1,1], (1+aLenExtension))
    lineEnd =  np.vstack((charAvg[-1-indexLength:-1,:], np.vstack(np.asarray([[lineend_0 , lineend_1]]))))
        
    #%Find intercepts to divide line using Poly
    _, iIntStart = poly.polyxpoly(closedEnvelope, lineStart)
    iIntStart = np.floor(iIntStart[0,0]).astype(int)
      
    _, iIntEnd = poly.polyxpoly(closedEnvelope, lineEnd)
    iIntEnd = np.floor(iIntEnd[-1,0]).astype(int)

  #% To divide inner or outer corridors, first determine if polygon is clockwise
  #% or counter-clockwise. Then, based on which index is large, separate out
  #% inner and outer corridor based on which intercept index is larger. 
  if poly.ispolycw(envelope):
    if iIntStart > iIntEnd:
      outerCorr = np.vstack([envelope[iIntStart:,:],envelope[:iIntEnd,:]])
      innerCorr = envelope[iIntEnd:iIntStart,:]
    else:
      outerCorr = envelope[iIntStart:iIntEnd,:]
      innerCorr = np.vstack([envelope[iIntEnd:,:], envelope[:iIntStart,:]])
  else:
    if iIntStart > iIntEnd:
      innerCorr = np.vstack([envelope[iIntStart:,:], envelope[:iIntEnd,:]])
      outerCorr = envelope[iIntEnd:iIntStart,:]
    else:
      innerCorr = envelope[iIntStart:iIntEnd,:]
      outerCorr = np.vstack([envelope[iIntEnd:,:],envelope[:iIntStart,:]])

  #% Resample corridors. Use nResamplePoints. Because corridors are
  #% non-monotonic, arc-length method discussed above is used. 
  #% Start with inner corridor. Magnitudes are being normalized.
  segments = np.sqrt(((innerCorr[0:-1,0]-innerCorr[1:,0]) / np.max(innerCorr[:,0])) **2 + ((innerCorr[0:-1,1]-innerCorr[1:,1])/np.max(innerCorr[:,1])) **2)
  alen = np.cumsum(np.concatenate([[0],segments]))
  alenResamp = np.linspace(0,np.max(alen), num = nResamplePoints)
  alenResamp = np.transpose(alenResamp)
  innerCorr = np.column_stack([interpolate.interp1d(alen,innerCorr[:,0])(alenResamp), interpolate.interp1d(alen,innerCorr[:,1])(alenResamp)])

  #% Outer Corridor
  segments = np.sqrt(((outerCorr[0:-1,0]-outerCorr[1:,0]) / np.max(outerCorr[:,0])) **2 + ((outerCorr[0:-1,1]-outerCorr[1:,1])/np.max(outerCorr[:,1])) **2)
  alen = np.cumsum(np.concatenate([[0],segments]))
  alenResamp = np.linspace(0,np.max(alen), num = nResamplePoints)
  alenResamp = np.transpose(alenResamp)
  outerCorr = np.column_stack([interpolate.interp1d(alen,outerCorr[:,0])(alenResamp), interpolate.interp1d(alen,outerCorr[:,1])(alenResamp)])


  #%% Draw extension lines and sampling points to MS plot
  if (Diagnostics == 'detailed'):
    figd, axd = plt.subplots(1, 2, figsize = (12,4), dpi = 1200)
    axd[0].plot(charAvg[:,0], charAvg[:,1], label= 'Char Avg', c = 'black')
    ellipse_xy = list(zip(charAvg[:,0], charAvg[:,1]))
    for iPoint in range(nResamplePoints):
      ellipse = Ellipse(ellipse_xy[iPoint], stdevData[iPoint,0] * EllipseKFact*2, stdevData[iPoint,1] * EllipseKFact*2, angle = 0, edgecolor='grey', lw=0.6, facecolor='none')
      axd[0].add_artist(ellipse)
    for iSignal in inputSignals.keys():
      axd[0].plot(inputSignals[iSignal][:,0], inputSignals[iSignal][:,1], label = iSignal)
    axd[0].plot(closedEnvelope[:,0], closedEnvelope[:,1], c = 'darkgreen', linewidth=0.5)

    axd[0].title.set_text('Char Avg and Ellipses')
    axd[0].set(xlabel='x-data', ylabel='y-data')
    axd[0].legend(loc='lower right')


    axd[1].scatter(xx[:],yy[:], 0.25, zz[:]>=1, )
    axd[1].title.set_text('Corridor Extraction')
    axd[1].set(xlabel='x-data', ylabel='y-data')
    figd.savefig('outputs/DetailedDiagnostics.png')

  # print average and corridors to .csv file
  if resultsToFile:
    output = np.column_stack([charAvg,innerCorr,outerCorr])
    fmt = ",".join(["%s"] + ["%10.6e"] * (output.shape[1]-1))
    np.savetxt("outputs/ArcGen Output.csv", output, fmt=fmt, header='Average Corridor, , Inner Corridor, , Outer Corridor, , \n x-axis, y-axis, x-axis, y-axis, x-axis, y-axis', comments='')
    fig = plt.figure(figsize= (6,4), dpi=1200)
    if NormalizeSignals == 'on':
      for iSignal in inputSignals.keys():
        #% Resulting array is normalized arc - length, resampled x, resam.y
        plt.plot(normalizedSignal[iSignal][:,1], normalizedSignal[iSignal][:,2], label = 'Input Signals', c ='grey', lw=1)
      plt.title('ArcGen - Normalization')
    else:
      for iSignal in inputSignals.keys():

        plt.plot(inputSignals[iSignal][:,0], inputSignals[iSignal][:,1], label = 'Input Signals', c = 'grey', lw = 1)
      plt.title('ArcGen - No Normalization')

    plt.plot(charAvg[:,0], charAvg[:,1], label = 'Average - ARCGen', c ='black', lw = 1.5)
    plt.plot(innerCorr[:,0], innerCorr[:,1], label = 'Inner Corridors', c='gold', lw = 1.5)
    plt.plot(outerCorr[:,0], outerCorr[:,1], label = 'Outer Corridors', c ='goldenrod', lw = 1.5)
    handles, labels = plt.gca().get_legend_handles_labels()
    labels, ids = np.unique(labels, return_index=True)
    handles = [handles[i] for i in ids]
    fig.legend(handles, labels)
    plt.xlabel('x - axis')
    # Set the y axis label of the current axis.
    plt.ylabel('y - axis')
    # Set a title of the current axes.
    # show a legend on the plot
    # Display a figure.
    fig.savefig('outputs/ARCGen - Corridors and Signal Avg.png')

  # Create output of processed signals
  processedSignals = []
  for key in inputSignals.keys():
    if nWarpCtrlPts > 0:
      tempDict = {
        "resampledSignals": normalizedSignal[key],
        "warpControlPoints": warpControlPoints[key],
      }
    else:
      tempDict = {
        "resampledSignals": normalizedSignal[key],
      }
    processedSignals.append(tempDict)

  # Create debug dictionary
  debugData = {
    "charAvg": charAvg,
    "stDev": stdev,
    "preWarpCorrArray": preWarpCorrArray,
    "preWarpMeanCorr": preWarpCorrArray,
    "warpedCorrArray": warpedCorrArray,
    "warpedMeanCorrScore": warpedMeanCorrScore,
  }

  return charAvg, innerCorr, outerCorr, processedSignals, debugData
	
## External Functions
#%% Function used to evaluate correlation score between signals
def evalCorrScore(signalsX, signalsY):
  #% Correlation score taken from the work of Nusholtz et al. (2009)
  #% Compute cross-correlation matrix of all signals to each other
  corrMatX = np.corrcoef(signalsX, rowvar=False)
  corrMatY = np.corrcoef(signalsY, rowvar=False)
  
  #% Convert matrices to a single score
  nSignal = len(corrMatX)
  corrScoreX = (1/(nSignal*(nSignal-1)))*(np.sum(np.sum(corrMatX))-nSignal)
  corrScoreY = (1/(nSignal*(nSignal-1)))*(np.sum(np.sum(corrMatY))-nSignal)
  
  #% Compute a single metric for optimization purposes. Using simple mean
  meanCorrScore = 0.5*(corrScoreX+corrScoreY)
  corrScoreArray = np.column_stack((corrScoreX, corrScoreY))
  return meanCorrScore, corrScoreArray


#%% Function used to compute objective for optimization
def warpingObjective(optimWarp,nCtrlPts,inputSignals,WarpingPenalty,nResamplePoints):
  #% Control points are equally spaced in arc-length.
  #% optimwarp is a column vector with first warped control point in the
  #% first nSignal indices, then 2nd control point in the next nSignal indices
  nSignal = len(inputSignals)
  warpArray = optimWarp.reshape((2*nSignal, nCtrlPts), order='F')
  
  #% Compute a warping penalty
  penaltyScore = warpingPenalty(warpArray, WarpingPenalty, nResamplePoints, nSignal)
  penaltyScore = np.mean(penaltyScore, axis=0)

  # % Perform warping
  _, signalsX, signalsY = warpArcLength(warpArray, inputSignals, nResamplePoints)
  
  #% Compute correlation score
  corrScore, _ = evalCorrScore(signalsX, signalsY)
  
  #% corrScore is a maximization goal. Turn into a minimization goal
  optScore = 1-corrScore+penaltyScore
  return optScore


#%% Function used to warp arc-length function
def warpArcLength(warpArray, inputSignals, nResamplePoints):

  #% Warp array: each row is warping points for an input signal, each column
  #% is warped point. Control points are interpolated  on [0,1] assuming
  #% equal spacing.
  nSignal = len(inputSignals)

  #% Initialize matrices
  signalsX = np.zeros((nResamplePoints, nSignal))
  signalsY = np.zeros((nResamplePoints, nSignal))
  warpedSignals = {}
  
  for iSignal, keys in enumerate(inputSignals):
    signal = inputSignals[keys]
    lmCtrlPts = np.concatenate([[0], warpArray[iSignal + nSignal, :], [1]])

    #% prepend 0 and append 1 to warp points for this signal to create valid
    #% control points.
    warpedCtrlPts = np.concatenate([[0], warpArray[iSignal,:], [1]])

    #% Construct warping function using SLM. This warps lmAlen to shiftAlen.
    #% Use warping fuction to map computed arc-lengths onto the shifted
    #% system. use built-in pchip function. This is a peicewise monotonic
    #% cubic spline. Signifincantly faster than SLM.
    warpedNormAlen = interpolate.pchip(lmCtrlPts,warpedCtrlPts)(signal[:,3])

    #% Now uniformly resample normalzied arc-length
    resamNormwarpedAlen = np.linspace(0, 1, num = nResamplePoints)
    resampX = interpolate.interp1d(warpedNormAlen, signal[:, 0], kind='linear', fill_value="extrapolate")(resamNormwarpedAlen)
    resampY = interpolate.interp1d(warpedNormAlen, signal[:, 1], kind='linear', fill_value="extrapolate")(resamNormwarpedAlen)

    #% Assign to array for correlation calc
    signalsX[:, iSignal] = resampX
    signalsY[:, iSignal] = resampY

    #% Assemble a cell array containing arrays of resampled signals. Similar
    #% to 'normalizedSignal' in 'inputSignals' structure
    warpedSignals[keys] = np.column_stack((resamNormwarpedAlen, resampX, resampY))
  
  return warpedSignals, signalsX, signalsY


#%% Penalty function to prevent plateaus and extreme divergence in warping functions
def warpingPenalty(warpArray, WarpingPenalty, nResamplePoints, nSignal):
#% Compute an array of penalty scores based on MSE between linear, unwarped
#% arc-length and warped arc-length. Aim is to help prevent plateauing.
  nSignals = nSignal

  penaltyScores = np.zeros((nSignals));
  unwarpedAlen = np.linspace(0, 1, num = nResamplePoints);
  
  for iSignal in range(nSignals):
    interpX = np.concatenate([[0],warpArray[iSignal+nSignals,:],[1]], axis=None, dtype='float')
    interpY = np.concatenate([[0],warpArray[iSignal,:],[1]], axis=None, dtype='float')
    interpResults = interpolate.pchip(interpX, interpY, axis=0)(unwarpedAlen)
    penaltyScores[iSignal] = np.sum(((unwarpedAlen - interpResults)**2))
  
  penaltyScores = WarpingPenalty * penaltyScores

  return penaltyScores


#%% helper function to perform linear interpolation to an isovalue of 1 only
def interpIso(x1, y1, x2, y2):
  val = x1+(x2-x1)*(1-y1)/(y2-y1)  
  return val

def interpLin(x1, y1, x2,  y2, xq):
  return y1 + (xq-x1)*(y2-y1)/(x2-x1)


#%% returns  indices and values of nonzero elements
def find(array, string = None):
  #making sure that a np array is used
  array = np.array(array)
  n = array.ndim
  if n ==0:
    return []
  elif n == 1:
    result = np.argwhere(array)
  elif n == 2:
    array1d = array.ravel(order='F')
    result = np.argwhere(array1d)
  else:
    raise ValueError("Array dimensions not supported")

  if string == 'first' and result.size != 0:
    return np.sort(result, axis = None)[0]

  elif string == 'last' and result.size != 0:
    return np.sort(result, axis = None)[-1]

  else:
    return np.array(result).astype(int)
