import arcgen.polygonFunctions as poly
from arcgen.uniquetol import uniquetol
import numpy as np
from numpy import inf
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy import interpolate
from scipy import optimize
import os
import warnings
warnings.filterwarnings("ignore")

def arcgen(inputData,
           nResamplePoints = 250,
           CorridorRes = 250,
           Diagnostics = 'off',
           NormalizeSignals = 'on',
           EllipseKFact = 1,
           MinCorridorWidth = 0,
           nWarpCtrlPts = 0,
           WarpingPenalty = 1e-2,
           resultsToFile = False):
    """
    ARCGen: Calculation of a characteristic average and statistical response
    corridors using arc-length re-parameterization and signal registration. 
    ARCGen can accomandate a very wide range of input signals, including non-
    monotonic nad hystertic signals.

    If you intend to publish work using ARCGen, please use the following 
    citation:
        Hartlen, D.C. & Cronin, D.S. (2022) "Arc-Length Re-Parametrization and 
            Signal Registration to Determine a Characteristic Average and 
            Statistical Response Corridors of Biomechanical Data." Frontiers in
             Bioengineering and Biotechnology 10:843148. 
             doi: 10.3389/fbioe.2022.843148

    ARCGen is released under a GNU GPL v3 license. No warranty or support is
    provided. The authors hold no responsibility for the validity, accuracy,
    or applicability of any results obtained from this code.

    Parameters
    ----------
    inputData: 
        One of three input formats:
        str: a path to a csv file containing all two-column signals 
        list of np.ndarrays: useful when signals are processed prior to ARCGen
        list of dictionaries: allows one to specify specimen ID. Dictionary
            format: {'data': np.ndarray, 'specId', str}
    nResamplePoints: int
        Number of points used to resample input signals during arc-length
            re-parameterization. Default: 250
    CorridorRes: int
        Size of grid (along one size) used for extracting corridors via a 
            marching squares algorithm. Grid extends 120% beyond the extreme
            corridor extents. Number of points in the grid is CorridorRes^2. 
            Default: 250
    nWarpCtrlPts: int
        Number of interior control points used during signal registration. A 
            value of 0 (default) disables signal registration
    warpingPenalty: float
        Penalty factor used to control warping during signal registration. 
            Default: 1e-2
    Diagnostics: str
        Controls if debugging plots are produced. Options are: 'off' (default),
            'on', 'detailed'
    NormalizeSignals: str
        Cotnrols if signal scaling is performed to prevent differences in 
            extent between ordinate and abscissca from monopolizing arc-length.
            Options: 'on' (Default and highly recommended), 'off'
    EllipseKFact: float
        Scales the size of elliptical confidence regions produced at each point
            of the characteristic average, hence scaling corridor extent. 
            Default: 1.0
    MinCorridorWidth: float
        Enforces minimum corridor width as a factor of maximum standard
            deviation. Any st. dev. less than MinCorridorWidth*max(st. dev.) is 
            replaced with MinCorridorWidth*max(st. dev.). Ordinate and 
            abscissca are handled independently. Default: 1.0.
    resultsToFile: bool
        If True, arcgen() will write characteristic average and corridors to a 
            csv file as well as producing a plot of input signals, average and
            corridors. These are output to the folder 'outputs' in the current 
            working directory. Default: False

    Returns
    -------
    charAvg: np.ndarray
        Computed characteristic average (nResamplePoints x 2 array)
    innCorr: np.ndarray
        Computed inner corridor (nResamplePoints x 2 array)
    outCorr: np.ndarray
        Computed outer corridor (nResamplePoints x 2 array)
    processedSignals: list of dictionaries
        List with as many entries each containing a dictionary with data and
            information used when processing input signals. Dictionary format
            is {'data': np.ndarray, 'specId': str, 'xMax': float,
            'xMin': float, 'yMax': float, 'yMin': float, 'maxAlen': float, 
            'resampled': np.ndarray, 'warpControlPoints': np.ndarray}
    debug: dictionary
        Contains information useful for debugging purposes. Dictionary format:
            {'charAvg': np.ndarray, 'stDev': np.ndarray, 
            'preWarpCorrArray': np.ndarray, 'preWarpMeanCorr': float, 
            'warpedCorrArray': np.ndarray, 'warpedMeanCorrScore': float}
    
    Notes
    -----
    Please refer to package readme, repository 
    (https://github.com/IMMC-UWaterloo/ARCGen-Python), or provided article for
    a detailed description of how ARCGen functions

    It is common to see issues when running this function if the number of
    resampling points or corridor resolution is too sparse or signals
    exhibit significant variablity not accounted for through signal
    registration. This tends to manifest in either truncated corridors or the
    code termininating in an error. Often increasing resampling points or
    corridor resolution. Turning 'Diagnostics' to 'detailed' can help
    identify these issues.

    Computed corridors will often not extend all the way to the shared origin
    of input signals. This is because small low st. dev. at this shared point
    is too low to be captured during corridor extraction with the marching
    squares algorithm. There are two solutions to this problem. First, one
    could force a minimum corridors size using the 'MinCorridorWidth' option.
    Second, one could manually extend corridors in post-processing.

    The python version of ARCGen is generally updated less frequently than the
    original MATLAB version of ARCGen. Active developments are implimented
    in MATLAB prior to being ported to python. The MATLAB version of ARCGen 
    can be found at https://github.com/IMMC-UWaterloo/ARCGen

    Very special thanks is due to Ahmed Ibrahim, my fellow Ph.D. candidate
    at UWaterloo, for his help in getting this port started. 

    Copyright (c) 2022 Devon C. Hartlen
    """

    # Creating a directory for results
    if (Diagnostics == 'on') or (Diagnostics == 'detailed') or (resultsToFile):
        if not os.path.exists('outputs'):
            os.makedirs('outputs')

    # Create 'inputSignals' as a variable of convenience to keep track of all
    # pertenant data. Will consist of a list of dictionaries. 
    inputSignals = []

    # Handle input data for various cases: Path, standard list, list of dicts
    # First case: if inputData is a path
    if isinstance(inputData, str):
        # Load from csv. empty values in shorter signals are given as 'nan'
        dataframe = np.genfromtxt(inputData, delimiter=',', encoding=None)
        numberRows, numberCols = dataframe.shape
        # Error check
        if numberCols % 2 == 0:
            nSignal = int(numberCols/2)
        else:
            raise ValueError("The number of columns in the csv" + \
                " file is not even")    

        # Store input signals as list of arrays
        for iSig in range(nSignal):
            temp = dataframe[:, (2*iSig):(2*iSig+2)]
            # Remove 'nan' entries
            indexNan = np.logical_not(np.isnan(temp[:,0]))
            # add to dictionary
            inputSignals.append({'data': temp[indexNan,:], 
                'specId': 'Signal '+str(iSig+1) })
    
    # Now check for the other two input cases. Both start with a list check
    elif isinstance(inputData, list):
        # Case 2: list of np.ndarrays
        if isinstance(inputData[0], np.ndarray):
            nSignal = len(inputData)
            for iSig in range(len(inputData)):
                inputSignals.append({'data': inputData[iSig], 
                    'specId': 'Signal '+str(iSig+1)})
        
        # Case 3: list of dictionaries
        elif isinstance(inputData[0], dict):
            nSignal = len(inputData)
            for entry in inputData:
                inputSignals.append({'data': entry['data'], 
                    'specId': entry['specId']})
        else:
            raise ValueError('Incorrectly formatted input list')
    else:
        raise ValueError('Input is not a path or or correctly formatted')

    # Compute arclength based on input signal datapoints
    # Do not perform normalization
    if NormalizeSignals == 'off':
        for signal in inputSignals:
            temp = signal['data'].copy()  # Temporary for convenience

            # Compute arc - length between each data point
            segments = np.sqrt((temp[0:-1,0] - temp[1:,0])**2 
                                + (temp[0:-1, 1] - temp[1:, 1])**2)
            segments = np.append([0], segments)

            # Append cumulative arc length to data array
            alen = np.cumsum(segments)

            # Compute normalized arc - length
            signal['maxAlen'] = alen[-1]
            signal['data'] = np.column_stack((signal['data'], alen))
            signal['data'] = np.column_stack((signal['data'], 
                                              alen/signal['maxAlen']))

            # Determine max[x, y] data
            signal['xMax'], signal['yMax'] = temp.max(axis=0)

            # Remove spurious duplicates based on arc-length
            _,uniqueIndex = np.unique(signal['data'][:,3], axis=0, 
                                      return_index=True)
            signal['data'] = signal['data'][uniqueIndex,:]

    # Perform magnitude scaling based on bounding box
    if NormalizeSignals =='on':
        # Determine bounding box of individual signals
        for signal in inputSignals:
            temp = signal['data'].copy()  # Temporary for convenience
            # Determine min[x, y] and max[x, y]data
            signal['xMin'], signal['yMin'] = temp.min(axis=0)
            signal['xMax'], signal['yMax'] = temp.max(axis=0)
        # Compute average bounding box. Using shorthand to access list of dicts 
        xBound = [sum([signal['xMin'] for signal in inputSignals])/nSignal, 
            sum([signal['xMax'] for signal in inputSignals])/nSignal]
        yBound = [sum([signal['yMin'] for signal in inputSignals])/nSignal, 
            sum([signal['yMax'] for signal in inputSignals])/nSignal]

        # Normalize the axis of each signal, then do arc-length calcs
        for signal in inputSignals:
            # This needs to be .copy(), otherwise scaling effects inputSignals 
            temp = signal['data'].copy() # Temporary for convenience

            # Normalize from bounding box to [-1,1]
            temp[:,0] = temp[:,0] / (xBound[1]-xBound[0])
            temp[:,1] = temp[:,1] / (yBound[1]-yBound[0])

            # Compute arc - length between each data point
            segments = np.sqrt((temp[0:-1, 0] - temp[1:, 0]) ** 2 
                                + (temp[0:-1, 1] - temp[1:, 1])** 2)
            segments = np.append([0], segments)

            # Append cumulative arc length to data array
            alen = np.cumsum(segments)

            # Compute normalized arc - length
            signal['maxAlen'] = alen[-1]
            signal['data'] = np.column_stack((signal['data'], alen))
            signal['data'] = np.column_stack((signal['data'], 
                                              alen/signal['maxAlen']))

            # Remove spurious duplicates based on arc-length
            _,uniqueIndex = np.unique(signal['data'][:,3], axis=0,
                                      return_index=True)
            signal['data'] = signal['data'][uniqueIndex,:]

    # Uniformly resample x & y w.r.t. arc-length, signal['data'][:,3]
    for signal in inputSignals:
        normAlen = np.linspace(0, signal['data'][-1, 3], num = nResamplePoints)

        resampX = interpolate.interp1d(signal['data'][:, 3],
            signal['data'][:, 0])(normAlen)
        resampY = interpolate.interp1d(signal['data'][:, 3], 
            signal['data'][:, 1])(normAlen)

        # Resulting array is normalized arc - length, resampled x, resam.y
        signal['resampled'] = np.column_stack((normAlen, resampX, resampY))

    # For each resampled point, calc average and st. dev.
    # Initialize arrays
    charAvg = np.zeros((nResamplePoints, 2))
    stdevData = np.zeros((nResamplePoints, 2))
    for iPoints in range(nResamplePoints):
        temp = np.empty([0, 2])
        for signal in inputSignals:
            temp = np.vstack((temp, signal['resampled'][iPoints,1:3]))
        charAvg[iPoints,:]= temp.mean(axis=0)
        stdevData[iPoints,:] = temp.std(axis=0)

    # Align normalized arc-length signals based on minimized correlation.
    # Enabled by option 'nWarpCtrlPts'. If 0, skip alignment.
    if nWarpCtrlPts > 0:
        # Assemble signal matrices prior to correlation
        signalX = np.zeros((nResamplePoints, nSignal))
        signalY = np.zeros((nResamplePoints, nSignal))
        for iSig, signal in enumerate(inputSignals):
            signalX[:,iSig] = signal['resampled'][:, 1]
            signalY[:,iSig] = signal['resampled'][:, 2]
        meanCorrScore, corrArray = evalCorrScore(signalX,signalY)

        # Assign pre-optimized correlation scores to debug structure
        preWarpCorrArray = corrArray
        preWarpMeanCorrScore = meanCorrScore

        # Optimize warp points for arbitrary n warping points. Build bounds,
        # constraints, and x0s
        nWarp = nWarpCtrlPts
        # nWarp == 1 is a special case as inequalites aren't needed
        if nWarp == 1:   
            x0 = 0.5 * np.ones((nSignal*2))
            lb = 0.15 * np.ones((nSignal*2))
            ub = 0.85 * np.ones((nSignal*2))
            A = []
            b = []
            linConstraints = []
        elif nWarp >= 15:
            raise ValueError('Specifying more than 15 interior warping points'\
                ' is not supported')
        else:
            x0 = np.zeros((nWarp * (nSignal * 2)))

            for i in range(nWarp):
                x0[np.array([(i * nSignal) + np.arange(0, nSignal)
                    + (i * (nSignal))])] \
                    = ((i+1)/(nWarp+1)) * np.ones((nSignal))
                x0[np.array([(i * nSignal) + np.arange(0, nSignal)
                    + ((i+1) * nSignal)])] \
                    = ((i+1)/(nWarp+1)) * np.ones((nSignal))

            lb = 0.05 * np.ones((nWarp*(nSignal*2)))
            ub = 0.95 * np.ones((nWarp*(nSignal*2)))
            A = np.zeros([((nWarp-1) * (nSignal * 2)), 
                           (nWarp * (nSignal * 2))]) 
             # Force some separation between warped points
            b = -0.05 * np.ones(((nWarp-1) * (nSignal * 2))) 
            
            for iSignal in range(nSignal * 2):
                for iWarp in range(nWarp-1):
                    A[iSignal + (iWarp) * (nSignal * 2),
                        iSignal + ((iWarp) * (nSignal * 2))] = 1
                    A[iSignal + (iWarp) * (nSignal * 2),
                        iSignal + ((iWarp+1) * (nSignal * 2))] = -1
            
            linConstraints = optimize.LinearConstraint(A,
                -np.inf * np.ones(((nWarp-1) * (nSignal * 2))), b)
        
        # Execute optimization and compute warped signals
        bounds = list(zip(lb, ub))  
        optWarpArrayX = \
            optimize.minimize(warpingObjective, 
                              x0, 
                              args=(nWarp,
                                    inputSignals, 
                                    WarpingPenalty,
                                    nResamplePoints),
                              method = 'SLSQP',
                              bounds=bounds, 
                              constraints=linConstraints,
                              options={ 'disp': False})
        optWarpArray = optWarpArrayX.x
        optWarpArray = optWarpArray.reshape((nSignal*2,nWarp), order='F')

        warpedSignals, signalX, signalY = warpArcLength(optWarpArray,
                                                        inputSignals,
                                                        nResamplePoints)

        # Compute correlation score
        meanCorrScore, corrArray = evalCorrScore(signalX, signalY)
        # Assign warped correlation scores to debug structure
        warpedCorrArray = corrArray
        warpedMeanCorrScore = meanCorrScore
                
        # Replace 'resampled' in 'inputSignals' and compute average and 
        # standard deviation.
        for iSig, signal in enumerate(inputSignals):
            signal['resampled'] = warpedSignals[iSig]
            tempRow1 = np.concatenate([[0],optWarpArray[iSig+nSignal,:],[1]])
            tempRow2 = np.concatenate([[0],optWarpArray[iSig,:],[1]])
            signal['warpControlPoints'] = \
                np.column_stack([tempRow1.T, tempRow2.T])
        
        # Recalculate average and standard devation based on warped arc-length
        for iPt in range(nResamplePoints):
            temp = np.empty((nSignal, 2)) # probably cleaner way to do this.
            for iSig, signal in enumerate(inputSignals):
                # collect specified point from each signal
                temp[iSig,:] = signal['resampled'][iPt,1:3]
            charAvg[iPt,:] = np.mean(temp, axis = 0)
            stdevData[iPt,:] = np.std(temp, axis = 0,  ddof=1)

    # Clamp minimum corridor width. Disabled if 'MinCorridorWidth' == 0
    # Include influence of corridor scaling factor 'EllipseKFact'
    if MinCorridorWidth > 0:
        # Replace any stDevData below maximum st.dev. * 'MinCorridorWidth'
        stdevData[:,0] = np.where(stdevData[:,0] < 
            (MinCorridorWidth * np.max(stdevData[:,0]) * EllipseKFact), 
            (MinCorridorWidth * np.max(stdevData[:,0]) * EllipseKFact),
            stdevData[:,0])
        stdevData[:,1] = np.where(stdevData[:,1] < 
            (MinCorridorWidth * np.max(stdevData[:,1]) * EllipseKFact),
            (MinCorridorWidth * np.max(stdevData[:,1]) * EllipseKFact),
            stdevData[:,1])

    # Diagnostic: Plot normalized signals and St. Devs. 
    if Diagnostics == 'on' or Diagnostics == 'detailed':
        # Plot normalized x,y 
        fig, ax = plt.subplots(2, 2, figsize = (8,6), dpi = 300)
        fig.suptitle('Arc-length Discretized Normalized Signals')
        ax[0][0].title.set_text('Arc-length Discretized Normalized Signals')
        ax[0][1].title.set_text('Warping Functions')
        ax[1][0].title.set_text('Average and St.Dev. of X-Data')
        ax[1][1].title.set_text('Average and St.Dev. of Y-Data')

        for iSig, signal in enumerate(inputSignals):
            # Plot resampled signalls
            ax[0][0].plot(signal['resampled'][:,1], signal['resampled'][:,2], 
                          label= signal['specId'])

            if nWarpCtrlPts > 0:
                interpResults = \
                    interpolate.pchip(signal['warpControlPoints'][:,0], 
                                      signal['warpControlPoints'][:,1], 
                                      axis=0)(signal['data'][:,3])

                ax[0][1].plot(signal['data'][:,3], interpResults, 
                              label = signal['specId'])
                ax[0][1].plot(signal['warpControlPoints'][:,0], 
                              signal['warpControlPoints'][:,1],
                              lw=0.0, marker='x', markersize=6)
            else: 
                ax[0][1].title.set_text('No Warping Performed')

            ax[1][0].plot(signal['resampled'][:,0], signal['resampled'][:,1], 
                          label= signal['specId'])
            ax[1][1].plot(signal['resampled'][:,0], signal['resampled'][:,2], 
                          label= signal['specId'])

        ax[1][0].errorbar(inputSignals[0]['resampled'][:,0], charAvg[:,0], 
                          yerr=stdevData[:,0], color='darkgrey')
        ax[1][1].errorbar(inputSignals[0]['resampled'][:,0], charAvg[:,1], 
                          yerr=stdevData[:,1], color='darkgrey')
            
        ax[0][1].plot([0,1],[0,1])
        ax[0][0].set(xlabel='x-data', ylabel='y-data')
        ax[0][1].set(xlabel='Unwarped Normalized Arc-length',
                     ylabel='Warped Normalized Arc-length')
        ax[1][0].set(xlabel='Normalized Arc-length', ylabel='x-data')
        ax[1][1].set(xlabel='Normalized Arc-length', ylabel='y-data')

        fig.tight_layout()
        fig.subplots_adjust(top=0.91)
        ax[0][0].legend(loc='lower right')
        fig.savefig('outputs/ARCGen Diagnostics.png')

    # Begin marching squares algorithm
    # Create grids based on upper and lower of characteristic average plus 120%
    # of maximum standard deviation
    scaleFact = 1.25 * EllipseKFact
    xx,yy = np.meshgrid(
        np.linspace((np.min(charAvg[:,0]) - scaleFact*np.max(stdevData[:,0])), 
                    (np.max(charAvg[:,0]) + scaleFact*np.max(stdevData[:,0])), 
                    num = CorridorRes),
        np.linspace((np.min(charAvg[:,1]) - scaleFact*np.max(stdevData[:,1])), 
                    (np.max(charAvg[:,1]) + scaleFact*np.max(stdevData[:,1])), 
                    num = CorridorRes),
        indexing='xy')
    # Using separate function here for the potential of future optimization
    zz = evaluateGrid(xx, yy, charAvg, stdevData, EllipseKFact)
            
    # The following segments is the marching squares algorith. The goal of this
    # algorithm is to find the zz=1 isoline, as this represents the outer
    # boundary of all elllipses. 
    #
    # Described in brief, this algorithm goes through each point, looking at
    # its and its neighbours values. There are only 16 configurations of these
    # squares or cells. Based on the configuration, add the appropriate line
    # segments. This method uses linear interpolation to increase accuracy. 
    # Initalize line segments for speed. This line may cause issues, as it
    # assumes maximum size. Bump up 10 if it does. 
    lineSegments = np.zeros((10 * max(nResamplePoints,CorridorRes), 4))

    iSeg = -1

    for iPt in range(CorridorRes-1):  # Rows (y-axis)
        for jPt in range(CorridorRes-1):   # Columns (x-axis)
            # Cell value definition
            #  1 -- 2 
            #  |    |
            #  |    |
            #  8 -- 4
            #
            # REMEMBER!!!! 
            # array(i,j) = array(rows, columns,) = array(y,x)
            
            # By carefully defining cell values and definitions, we can use
            # binary to simplify logic though a integer based switch case

            cellValue = int((1 * (zz[iPt,jPt]>1)) + 
                            (2 * (zz[iPt+1,jPt]>1)) + 
                            (4 * (zz[iPt+1,jPt+1]>1)) + 
                            (8 * (zz[iPt,jPt+1]>1)) + 1)
            
            if cellValue == 1:
                # No Vertices
                pass

            elif  cellValue == 2:
                # South-West
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    interpIso(xx[iPt,jPt], zz[iPt,jPt],
                               xx[iPt,jPt+1], zz[iPt,jPt+1]),
                    yy[iPt,jPt],
                    xx[iPt,jPt], 
                    interpIso(yy[iPt,jPt],zz[iPt,jPt],
                               yy[iPt+1,jPt],zz[iPt+1,jPt])]

            elif  cellValue == 3:
                # West-North
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    xx[iPt+1,jPt],
                    interpIso(yy[iPt,jPt], zz[iPt,jPt], 
                              yy[iPt+1,jPt], zz[iPt+1,jPt]),
                    interpIso(xx[iPt+1,jPt], zz[iPt+1,jPt], 
                              xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),
                    yy[iPt+1,jPt]]
            
            elif  cellValue == 4:
                # North-South
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    interpIso(xx[iPt,jPt], zz[iPt,jPt], 
                              xx[iPt,jPt+1], zz[iPt,jPt+1]),
                    yy[iPt,jPt],
                    interpIso(xx[iPt+1,jPt], zz[iPt+1,jPt],
                              xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),
                    yy[iPt+1,jPt]]
                        
            elif  cellValue == 5:
                # North-East
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    interpIso(xx[iPt+1,jPt], zz[iPt+1,jPt], 
                              xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),
                    yy[iPt+1,jPt+1],
                    xx[iPt+1,jPt+1], 
                    interpIso(yy[iPt+1,jPt+1], zz[iPt+1,jPt+1], 
                              yy[iPt,jPt+1], zz[iPt,jPt+1])]
                        
            elif  cellValue == 6:  # Ambiguous 
                centerVal = np.mean(list([zz[iPt,jPt], zz[iPt+1,jPt], 
                                          zz[iPt+1,jPt+1], zz[iPt, jPt+1]]))
                if centerVal >= 1:
                    # West-North
                    iSeg = iSeg+1
                    lineSegments[iSeg,:] = [
                        xx[iPt+1,jPt],
                        interpIso(yy[iPt,jPt], zz[iPt,jPt], 
                                  yy[iPt+1,jPt], zz[iPt+1,jPt]),
                        interpIso(xx[iPt+1,jPt], zz[iPt+1,jPt], 
                                  xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),
                        yy[iPt+1,jPt]]
                                
                    # South-East
                    iSeg = iSeg+1
                    lineSegments[iSeg,:] = [
                        interpIso(xx[iPt,jPt+1], zz[iPt,jPt+1], 
                                  xx[iPt,jPt], zz[iPt,jPt]),
                        yy[iPt,jPt+1],
                        xx[iPt,jPt+1],
                        interpIso(yy[iPt,jPt+1], zz[iPt,jPt+1], 
                                  yy[iPt+1,jPt+1], zz[iPt+1,jPt+1])]
                else:
                    # South-West
                    iSeg = iSeg+1
                    lineSegments[iSeg,:] = [
                        interpIso(xx[iPt,jPt], zz[iPt,jPt], 
                                  xx[iPt,jPt+1], zz[iPt,jPt+1]), 
                        yy[iPt,jPt],
                        xx[iPt,jPt], 
                        interpIso(yy[iPt,jPt], zz[iPt,jPt], 
                                  yy[iPt+1,jPt], zz[iPt+1,jPt])]
                                
                    # North-East
                    iSeg = iSeg+1
                    lineSegments[iSeg,:] = [
                        interpIso(xx[iPt+1,jPt], zz[iPt+1,jPt], 
                                  xx[iPt+1,jPt+1],zz[iPt+1,jPt+1]),
                        yy[iPt+1,jPt+1],
                        xx[iPt+1,jPt+1],
                        interpIso(yy[iPt+1,jPt+1], zz[iPt+1,jPt+1], 
                                  yy[iPt,jPt+1], zz[iPt,jPt+1])]

            elif  cellValue == 7:
                # West-East
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    xx[iPt,jPt],
                    interpIso(yy[iPt,jPt], zz[iPt,jPt],
                              yy[iPt+1,jPt], zz[iPt+1,jPt]),
                    xx[iPt,jPt+1],
                    interpIso(yy[iPt,jPt+1], zz[iPt,jPt+1], 
                              yy[iPt+1,jPt+1], zz[iPt+1,jPt+1])]
            
            elif  cellValue == 8:
                # South - East
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    interpIso(xx[iPt,jPt+1], zz[iPt,jPt+1], 
                              xx[iPt,jPt], zz[iPt,jPt]),
                    yy[iPt,jPt+1],
                    xx[iPt,jPt+1],
                    interpIso(yy[iPt,jPt+1], zz[iPt,jPt+1], 
                              yy[iPt+1,jPt+1], zz[iPt+1,jPt+1])]
                                                                            
            elif  cellValue == 9:
                # South - East
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    interpIso(xx[iPt,jPt+1], zz[iPt,jPt+1],
                              xx[iPt,jPt], zz[iPt,jPt]),
                    yy[iPt,jPt+1],
                    xx[iPt,jPt+1],
                    interpIso(yy[iPt,jPt+1], zz[iPt,jPt+1],
                              yy[iPt+1,jPt+1], zz[iPt+1,jPt+1])]

            elif  cellValue == 10:
                # West-East
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    xx[iPt,jPt],
                    interpIso(yy[iPt,jPt], zz[iPt,jPt],
                              yy[iPt+1,jPt], zz[iPt+1,jPt]),
                    xx[iPt,jPt+1],
                    interpIso(yy[iPt,jPt+1], zz[iPt,jPt+1], 
                              yy[iPt+1,jPt+1],zz[iPt+1,jPt+1])]

            elif  cellValue == 11: # Ambiguous
                centerVal = np.mean(list([zz[iPt,jPt], zz[iPt+1,jPt], 
                                          zz[iPt+1,jPt+1], zz[iPt, jPt+1]]))
                if centerVal >= 1:
                    # South-West
                    iSeg = iSeg+1
                    lineSegments[iSeg,:] = [
                        interpIso(xx[iPt,jPt], zz[iPt,jPt],
                                  xx[iPt,jPt+1], zz[iPt,jPt+1]),
                        yy[iPt,jPt],
                        xx[iPt,jPt], 
                        interpIso(yy[iPt,jPt], zz[iPt,jPt], 
                                  yy[iPt+1,jPt], zz[iPt+1,jPt])]

                    # North-East
                    iSeg = iSeg+1
                    lineSegments[iSeg,:] = [
                        interpIso(xx[iPt+1,jPt], zz[iPt+1,jPt],
                                  xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),
                        yy[iPt+1,jPt+1],
                        xx[iPt+1,jPt+1],
                        interpIso(yy[iPt+1,jPt+1], zz[iPt+1,jPt+1], 
                                  yy[iPt,jPt+1], zz[iPt,jPt+1])]
                else:
                    # West-North
                    iSeg = iSeg+1
                    lineSegments[iSeg,:] = [
                        xx[iPt+1,jPt],
                        interpIso(yy[iPt,jPt], zz[iPt,jPt], 
                                  yy[iPt+1,jPt], zz[iPt+1,jPt]),
                        interpIso(xx[iPt+1,jPt], zz[iPt+1,jPt],
                                  xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),
                        yy[iPt+1,jPt]]
                    # South-East
                    iSeg = iSeg+1
                    lineSegments[iSeg,:] = [
                        interpIso(xx[iPt,jPt+1], zz[iPt,jPt+1],
                                  xx[iPt,jPt], zz[iPt,jPt]), 
                        yy[iPt,jPt+1],
                        xx[iPt,jPt+1],
                        interpIso(yy[iPt,jPt+1], zz[iPt,jPt+1],
                                  yy[iPt+1,jPt+1], zz[iPt+1,jPt+1])]

            elif  cellValue == 12:
                # North-East
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    interpIso(xx[iPt+1,jPt], zz[iPt+1,jPt],
                              xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),
                    yy[iPt+1,jPt+1],
                    xx[iPt+1,jPt+1],
                    interpIso(yy[iPt+1,jPt+1], zz[iPt+1,jPt+1], 
                              yy[iPt,jPt+1], zz[iPt,jPt+1])]
            
            elif  cellValue == 13:
                # North-South
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    interpIso(xx[iPt,jPt], zz[iPt,jPt],
                              xx[iPt,jPt+1], zz[iPt,jPt+1]),
                    yy[iPt,jPt],
                    interpIso(xx[iPt+1,jPt], zz[iPt+1,jPt], 
                              xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),
                    yy[iPt+1,jPt]]
            
            elif  cellValue == 14:
                # West-North
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    xx[iPt+1,jPt],
                    interpIso(yy[iPt,jPt], zz[iPt,jPt], 
                              yy[iPt+1,jPt], zz[iPt+1,jPt]),
                    interpIso(xx[iPt+1,jPt], zz[iPt+1,jPt], 
                              xx[iPt+1,jPt+1], zz[iPt+1,jPt+1]),
                    yy[iPt+1,jPt]]

            elif  cellValue == 15:
                # South-West
                iSeg = iSeg+1
                lineSegments[iSeg,:] = [
                    interpIso(xx[iPt,jPt], zz[iPt,jPt], 
                              xx[iPt,jPt+1], zz[iPt,jPt+1]),
                    yy[iPt,jPt],
                    xx[iPt,jPt], 
                    interpIso(yy[iPt,jPt], zz[iPt,jPt],
                              yy[iPt+1,jPt], zz[iPt+1,jPt])]

            elif  cellValue == 16:
                pass
                # No vertices

    lineSegments = lineSegments[0:iSeg+1,:]

    # Extract list of unique vertices from line segmens
    vertices = np.vstack([lineSegments[:,:2],lineSegments[:,2:4]])
    vertices = uniquetol(vertices, 1e-12)
    index = np.zeros(len(vertices), dtype=bool)

    # Create a vertex connectivity table. The 1e-12 value is here because
    # floats don't round well and == will not work. 
    vertConn = np.zeros((len(lineSegments),2)).astype(int)

    for i in range(len(vertConn)):
        index = np.all(np.abs(lineSegments[:,:2] 
                              - vertices[i,:]) < 1e-12, axis = 1)
        vertConn[index, 0] = i
        index = np.all(np.abs(lineSegments[:,2:4] 
                              - vertices[i,:]) < 1e-12, axis = 1)
        vertConn[index, 1] = i

    # Start line segments sorting and envelope extraction
    nEnvelopes = 0
    allEnvelopes = np.zeros((len(vertConn),2))

    for i in range(len(vertConn)-1):
        j = i + 1 # helper index

        # save vertex to find 
        vertToFind = vertConn[i,1].astype(int)

        # Find connecting node
        vertConnToFind_ = np.any((vertConn[j:,:].astype(int) 
                                     == vertToFind.astype(int)), axis = 1)

        foundShiftedInd = find(vertConnToFind_.astype(int), 'first')

        # If we have found an index 
        if foundShiftedInd == 0 and vertConn[j,0] == vertToFind:
            continue

        elif foundShiftedInd == 0 and vertConn[j,1] == vertToFind:
            vertConn[j, [0, 1]] = vertConn[j, [1, 0]]
            continue
        
        elif foundShiftedInd > 0:
            foundInd = foundShiftedInd + j

            # swap found vert conn row with j row
            vertConn[[j,foundInd],:] = vertConn[[foundInd,j],:]
            if vertConn[j,0] == vertToFind:
                continue
            else:
                vertConn[j,[0,1]] = vertConn[j, [1,0]]
            
        # If we did not find an index, we either may have an open envelope or
        # envelope may be convex and loops back on itself. 
        elif foundShiftedInd.size == 0:
            # Check to see if we can find the first vertex in envelope
            # appearing again (check for closure)
            vertConn_ = allEnvelopes[nEnvelopes,0].astype(int)
            vertToFind = vertConn[vertConn_,1] 
            foundShiftedInd = find(np.any(vertConn[j:,:] 
                                      == vertToFind, axis = 1), 'first')
            
            # If we do not find an index, it means this envelope is complete
            # and manifold
            if foundShiftedInd.size == 0:
            # Assign indices to finish current envelope, initialize next
                allEnvelopes[nEnvelopes,1] = i
                nEnvelopes = nEnvelopes + 1
                allEnvelopes[nEnvelopes, 0] = j
            else:
                raise ValueError('Literal Edge Case')

    allEnvelopes[nEnvelopes,1] = j

    # Find largest envelope
    envInds = np.argmax((allEnvelopes[:,1]-allEnvelopes[:,0])).astype(int)
    # Convert indices in evelopes to array of (x,y)
    envInds = allEnvelopes[envInds, :]
    envelope = vertices[vertConn[envInds[0].astype(int):
                            envInds[1].astype(int)+1,0] ,:]

    # Divide the envelope into corridors. 
    # To break the largest envelop into inner and outer corridors, we need to
    # account for several edge cases. First, we test to see if there are any
    # intercepts of the characteristic average and the largest envelope. 
    closedEnvelope = np.vstack((envelope, envelope[0,:]))

    _, indexIntercept = poly.polyxpoly(closedEnvelope,charAvg)

    indexIntercept = np.floor(indexIntercept).astype(int)

    # If we find two intercepts, then we have no problem
    if indexIntercept.shape[0] >= 2:
        # Sort to order charAvg interecepts (2nd column), smaller charAvg 
        # intercepts are closer to start. 
        indexSort = np.argsort(indexIntercept[:,1], axis=0)
        iIntStart = indexIntercept[indexSort[0],0]
        iIntEnd = indexIntercept[indexSort[-1],0]

    # If we find only one intercept, we need to determine if the intercept is a
    # the start or end of the envelope. Then we need to extend the opposite
    # side of the characteristic average to intercept the envelope. 
    elif indexIntercept.shape[0] == 1:
        # Compute extension 
        aLenInterval = 1/nResamplePoints
        indexLength = round(0.2*len(charAvg))
            
        aLenExtension = (np.abs(aLenInterval/(charAvg[0,:]-charAvg[1,:])) 
                            * 1.1 * np.max(stdevData))
        aLenExtension[aLenExtension == inf] = 0
        aLenExtension = max(aLenExtension)
        
        # If the single found point is inside the envelope, the found intercept
        # is at the end. Therefore extend the start
        if poly.inpolygon(envelope, charAvg[indexIntercept[0,1],:] ): 
            iIntEnd = indexIntercept[-1,0]
            linestart_0 = interpLin(0, charAvg[0,0], 
                                    aLenInterval, charAvg[1,0], 
                                    -aLenExtension)
            linestart_1 = interpLin(0, charAvg[0,1], 
                                    aLenInterval, charAvg[1,1], 
                                    -aLenExtension)

            lineStart =  np.vstack((
                np.vstack(np.asarray([[linestart_0 , linestart_1]])), 
                charAvg[0:indexLength,:]))

        #Find intercepts to divide line using Poly
            _, iIntStart = poly.polyxpoly(closedEnvelope, lineStart)
            iIntStart = np.floor(iIntStart[0,0]).astype(int)

        # If the single found point is outside the envelope, the found
        # intercept is the start
        else:
            iIntStart = indexIntercept[0,0]
            # charAvg = np.asarray(charAvg)
            lineend_0 = interpLin(1-aLenInterval, charAvg[-2,0], 
                                  1, charAvg[-1,0], 
                                  (1+aLenExtension))
            lineend_1 = interpLin(1-aLenInterval, charAvg[-2,1], 
                                  1, charAvg[-1,1], 
                                  (1+aLenExtension))

            lineEnd =  np.vstack((
                charAvg[-1-indexLength:-1,:], 
                np.vstack(np.asarray([[lineend_0 , lineend_1]]))))
            
            # Find intercepts to divide line using Poly
            _, iIntEnd = poly.polyxpoly(closedEnvelope, lineEnd)
            iIntEnd = np.floor(iIntEnd[-1,0]).astype(int)
            
    # If we find no intercepts, we need to extend both sides of characteristic
    # average to intercept the envelop.
    else:
        aLenInterval = 1/nResamplePoints
        indexLength = round(0.2*len(charAvg))
        
        aLenExtension = (np.abs(aLenInterval/(charAvg[0,:]-charAvg[1,:])) 
                         * 1.1 * np.max(stdevData))
        aLenExtension[aLenExtension == inf] = 0
        aLenExtension = max(aLenExtension)

        linestart_0 = interpLin(0, charAvg[0,0], 
                                aLenInterval, charAvg[1,0], 
                                -aLenExtension)
        linestart_1 = interpLin(0, charAvg[0,1], 
                                aLenInterval, charAvg[1,1], 
                                -aLenExtension)
        lineStart =  np.vstack((
            np.vstack(np.asarray([[linestart_0 , linestart_1]])),
            charAvg[0:indexLength,:]))
        
        lineend_0 = interpLin(1-aLenInterval, charAvg[-2,0], 
                              1, charAvg[-1,0], 
                              (1+aLenExtension))
        lineend_1 = interpLin(1-aLenInterval, charAvg[-2,1], 
                              1, charAvg[-1,1], 
                              (1+aLenExtension))
        lineEnd =  np.vstack((
            charAvg[-1-indexLength:-1,:], 
            np.vstack(np.asarray([[lineend_0 , lineend_1]]))))
                
        # Find intercepts to divide line using polyxpoly
        _, iIntStart = poly.polyxpoly(closedEnvelope, lineStart)
        iIntStart = np.floor(iIntStart[0,0]).astype(int)
            
        _, iIntEnd = poly.polyxpoly(closedEnvelope, lineEnd)
        iIntEnd = np.floor(iIntEnd[-1,0]).astype(int)

    # To divide inner or outer corridors, first determine if polygon is 
    # clockwise or counter-clockwise. Then, based on which index is large, 
    # separate out inner and outer corridor based on which intercept index is 
    # larger. 
    if poly.ispolycw(envelope):
        if iIntStart > iIntEnd:
            outerCorr = np.vstack([envelope[iIntStart:,:],
                                   envelope[:iIntEnd,:]])
            innerCorr = envelope[iIntEnd:iIntStart,:]
        else:
            outerCorr = envelope[iIntStart:iIntEnd,:]
            innerCorr = np.vstack([envelope[iIntEnd:,:], 
                                   envelope[:iIntStart,:]])
    else:
        if iIntStart > iIntEnd:
            innerCorr = np.vstack([envelope[iIntStart:,:], 
                                   envelope[:iIntEnd,:]])
            outerCorr = envelope[iIntEnd:iIntStart,:]
        else:
            innerCorr = envelope[iIntStart:iIntEnd,:]
            outerCorr = np.vstack([envelope[iIntEnd:,:],
                                   envelope[:iIntStart,:]])

    # Resample corridors. Use nResamplePoints. Because corridors are
    # non-monotonic, arc-length method discussed above is used. 
    # Start with inner corridor. Magnitudes are being normalized.
    segments = np.sqrt(((innerCorr[0:-1,0] - innerCorr[1:,0]) 
                        / np.max(innerCorr[:,0]))**2 
                        + ((innerCorr[0:-1,1] - innerCorr[1:,1])
                        / np.max(innerCorr[:,1]))**2)
    alen = np.cumsum(np.concatenate([[0],segments]))
    alenResamp = np.linspace(0,np.max(alen), num = nResamplePoints)
    alenResamp = np.transpose(alenResamp)
    innerCorr = np.column_stack([
        interpolate.interp1d(alen,innerCorr[:,0])(alenResamp), 
        interpolate.interp1d(alen,innerCorr[:,1])(alenResamp)])

    # Outer Corridor
    segments = np.sqrt(((outerCorr[0:-1,0] - outerCorr[1:,0]) 
                         / np.max(outerCorr[:,0]))**2 
                         + ((outerCorr[0:-1,1] - outerCorr[1:,1])
                         / np.max(outerCorr[:,1]))**2)
    alen = np.cumsum(np.concatenate([[0],segments]))
    alenResamp = np.linspace(0,np.max(alen), num = nResamplePoints)
    alenResamp = np.transpose(alenResamp)
    outerCorr = np.column_stack([
        interpolate.interp1d(alen,outerCorr[:,0])(alenResamp), 
        interpolate.interp1d(alen,outerCorr[:,1])(alenResamp)])


    # Draw extension lines and sampling points to MS plot
    if (Diagnostics == 'detailed'):
        figd, axd = plt.subplots(1, 2, figsize = (12,4), dpi = 1200)
        axd[0].plot(charAvg[:,0], charAvg[:,1], label= 'Char Avg', c = 'black')
        ellipse_xy = list(zip(charAvg[:,0], charAvg[:,1]))
        for iPoint in range(nResamplePoints):
            ellipse = Ellipse(ellipse_xy[iPoint],
                             stdevData[iPoint,0] * EllipseKFact*2, 
                             stdevData[iPoint,1] * EllipseKFact*2, 
                             angle = 0, 
                             edgecolor='grey', 
                             lw=0.6, facecolor='none')
            axd[0].add_artist(ellipse)
        for signal in inputSignals:
            axd[0].plot(signal['data'][:,0], 
                        signal['data'][:,1], label = iSignal)
        axd[0].plot(closedEnvelope[:,0], closedEnvelope[:,1], 
                    c = 'darkgreen', linewidth=0.5)

        axd[0].title.set_text('Char Avg and Ellipses')
        axd[0].set(xlabel='x-data', ylabel='y-data')
        axd[0].legend(loc='lower right')


        axd[1].scatter(xx[:],yy[:], 0.25, zz[:]>=1, )
        axd[1].title.set_text('Corridor Extraction')
        axd[1].set(xlabel='x-data', ylabel='y-data')
        figd.savefig('outputs/Detailed Diagnostics.png')

    # print average and corridors to .csv file
    if resultsToFile:
        output = np.column_stack([charAvg,innerCorr,outerCorr])
        fmt = ",".join(["%s"] + ["%10.6e"] * (output.shape[1]-1))
        np.savetxt("outputs/ARCGen Output.csv", output, fmt=fmt, 
            header='Average Corridor, , Inner Corridor, , Outer Corridor, ,'
                    +'\n x-axis, y-axis, x-axis, y-axis, x-axis, y-axis', 
            comments='')

        fig = plt.figure(figsize= (6,4), dpi=1200)
        if NormalizeSignals == 'on':
            for signal in inputSignals:
                # Resulting array is normalized arc-length, resamp x, resam y
                plt.plot(signal['resampled'][:,1],
                         signal['resampled'][:,2], 
                         label = 'Input Signals', c ='grey', lw=1)
            plt.title('ArcGen - Normalization')
        else:
            for signal in inputSignals:
                plt.plot(signal['data'][:,0], 
                         signal['data'][:,1], 
                         label = 'Input Signals', c = 'grey', lw = 1)
            plt.title('ArcGen - No Normalization')

        plt.plot(charAvg[:,0], charAvg[:,1], 
                 label = 'Average - ARCGen', c ='black', lw = 1.5)
        plt.plot(innerCorr[:,0], innerCorr[:,1], 
                 label = 'Inner Corridors', c='gold', lw = 1.5)
        plt.plot(outerCorr[:,0], outerCorr[:,1], 
                 label = 'Outer Corridors', c ='goldenrod', lw = 1.5)
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
    processedSignals = inputSignals

    # Create debug dictionary
    if nWarpCtrlPts > 0:
        debugData = {
            "charAvg": charAvg,
            "stDev": stdevData,
            "preWarpCorrArray": preWarpCorrArray,
            "preWarpMeanCorr": preWarpCorrArray,
            "warpedCorrArray": warpedCorrArray,
            "warpedMeanCorrScore": warpedMeanCorrScore,
        }
    else:
        debugData = {
            "charAvg": charAvg,
            "stDev": stdevData,
        }

    return charAvg, innerCorr, outerCorr, processedSignals, debugData
    

def evalCorrScore(signalsX, signalsY):
    """
    Compute the correlation score of signals for use in signal registration. 
    Correlation score used here is taken from the work of Nusholtz et al.
    (2009)

    Parameters:
    -----------
    signalsX: np.ndarray
        X-component of re-parameterized signals, each signal is a column
    signalsY: np.ndarray
        Y-component of re-parameterized signals, each signal is a column

    Returns:
    --------
    meanCorrScore: float
        Average of x and y correlation scores
    corrScoreArray: np.ndarray
        Array of x and y axis correlation scores

    Copyright (c) 2022 Devon C. Hartlen
    """

    # Compute cross-correlation matrix of all signals to each other
    corrMatX = np.corrcoef(signalsX, rowvar=False)
    corrMatY = np.corrcoef(signalsY, rowvar=False)
    
    # Convert matrices to a single score
    nSignal = len(corrMatX)
    corrScoreX = (1/(nSignal*(nSignal-1)))*(np.sum(np.sum(corrMatX))-nSignal)
    corrScoreY = (1/(nSignal*(nSignal-1)))*(np.sum(np.sum(corrMatY))-nSignal)
    
    # Compute a single metric for optimization purposes. Using simple mean
    meanCorrScore = 0.5*(corrScoreX+corrScoreY)
    corrScoreArray = np.column_stack((corrScoreX, corrScoreY))
    return meanCorrScore, corrScoreArray


def warpingObjective(optimWarp,nCtrlPts,inputSignals,
                     penaltyFactor,nResamplePoints):
    """
    Objective function used by optimization algorithm during signal 
    registration with the goal of computing optimal control point locations. 
    Does not perform warping, only evaluates objective for optimization

    Parameters
    ----------
    optimWarp: np.ndarray
        Array of warping control points for which the objective function is 
        computed for
    nCtrlPts: int
        Number of interior control points
    inputSignals: list of dictionaries
        The internal list of dictionaries used to keep track of signal
        information inside arcgen()
    penaltyFactor: float
        Penalty factor used to control warping
    nResamplePoints: int
        Number of points that re-parameterized signals contains.
    Returns
    -------
    optScore: float
        objective score used for optimization
    
    Copyright (c) 2022 Devon C. Hartlen
    """

    # optimwarp is a column vector with first warped control point in the
    # first nSignal indices, then 2nd control point in the next nSignal indices
    nSignal = len(inputSignals)
    warpArray = optimWarp.reshape((2*nSignal, nCtrlPts), order='F')
    
    # Compute a warping penalty
    penaltyScore = warpingPenalty(warpArray, penaltyFactor, 
                                  nResamplePoints, nSignal)
    penaltyScore = np.mean(penaltyScore, axis=0)

    # Perform warping
    _, signalsX, signalsY = warpArcLength(warpArray, 
                                          inputSignals, 
                                          nResamplePoints)
    
    # Compute correlation score
    corrScore, _ = evalCorrScore(signalsX, signalsY)
    
    # corrScore is a maximization goal. Turn into a minimization goal
    optScore = 1-corrScore+penaltyScore
    return optScore


def warpArcLength(warpArray, inputSignals, nResamplePoints):
    """
    Given an array of warping control points, perform arc-length warping on 
    re-parameterized signals.

    Parameters
    ----------
    warpArray: np.ndarray
        Array of interior warping control points for all signals
    inputSignals: np.ndarray
        The internal list of dictionaries used by arcgen() to keep track of 
        signals and signal information
    nResamplePoints: int
        Number of points that re-parameterized signals contain
    
    Returns
    -------
    warpedSignals: list of np.ndarrays
        Resampled signals warped w.r.t arc-length using the control points 
        specified in warpArray
    signalsX: np.ndarray
        X-component of warped signals, each signal is a column
    signalsY: np.ndarray
        Y-component of warped signals, each signal is a column

    Copyright (c) 2022 Devon C. Hartlen
    """

    # Warp array: each row is warping points for an input signal, each column
    # is warped point. Control points are interpolated  on [0,1] assuming
    # equal spacing.
    nSignal = len(inputSignals)

    # Initialize matrices
    signalsX = np.zeros((nResamplePoints, nSignal))
    signalsY = np.zeros((nResamplePoints, nSignal))
    # Initialize a list for warped signal (a list of np.arrays)
    warpedSignals = [None]*nSignal
    
    for iSig, signal in enumerate(inputSignals):
        temp = signal['data'].copy()

        # prepend 0 and append 1 to warp points for this signal to create valid
        # control points.
        lmCtrlPts = np.concatenate([[0], warpArray[iSig + nSignal, :], [1]])
        warpedCtrlPts = np.concatenate([[0], warpArray[iSig,:], [1]])

        # warping function based on monotonic peicewise cubic hermite splines
        warpedNormAlen = interpolate.pchip(lmCtrlPts, warpedCtrlPts)(temp[:,3])

        # Now uniformly resample normalzied arc-length
        resamNormwarpedAlen = np.linspace(0, 1, num = nResamplePoints)
        resampX = interpolate.interp1d(warpedNormAlen, temp[:, 0], 
                                       kind='linear', 
                                       fill_value="extrapolate")\
                                       (resamNormwarpedAlen)
        resampY = interpolate.interp1d(warpedNormAlen, temp[:, 1], 
                                       kind='linear', 
                                       fill_value="extrapolate")\
                                       (resamNormwarpedAlen)

        # Assign to array for correlation calc
        signalsX[:, iSig] = resampX
        signalsY[:, iSig] = resampY

        # Assemble a cell array containing arrays of resampled signals. Similar
        # to 'normalizedSignal' in 'inputSignals' structure
        warpedSignals[iSig] = np.column_stack(
            (resamNormwarpedAlen, resampX, resampY))
    
    return warpedSignals, signalsX, signalsY


def warpingPenalty(warpArray, penaltyFactor, nResamplePoints, nSignal):
    """
    Computes the combined warping penalty for all signals given the provided
    warping points and inputted penalty factor

    Parameters
    ----------
    warpArray: np.ndarray
        Array of interior warping control points for all signals
    penaltyFactor: float
        Penalty factor used to control warping
    nResamplePoints: int
        Number of points that re-parameterized signals contain
    nSignal: int
        Number of input signals

    Returns
    -------
    penaltyScores: list, array-like
        Individual penalty scores for each warping function inputted. Same 
        length as nSignal
    
    Copyright (c) 2022 Devon C. Hartlen
    """
    # Compute an array of penalty scores based on MSE between linear, unwarped
    # arc-length and warped arc-length. Aim is to help prevent plateauing.
    penaltyScores = np.zeros((nSignal));
    unwarpedAlen = np.linspace(0, 1, num = nResamplePoints);
    
    for iSignal in range(nSignal):
        interpX = np.concatenate([[0], warpArray[iSignal+nSignal,:], [1]], 
                                 axis=None, dtype='float')
        interpY = np.concatenate([[0], warpArray[iSignal,:], [1]], 
                                 axis=None, dtype='float')
        interpResults = interpolate.pchip(interpX, interpY, axis=0)\
                                          (unwarpedAlen)
        penaltyScores[iSignal] = np.sum(((unwarpedAlen - interpResults)**2))
    
    penaltyScores = penaltyFactor * penaltyScores

    return penaltyScores


def interpIso(x1, y1, x2, y2):
    """
    Performs linear interpolation to find where x is equal to 1.0. Only used 
    for marching squares algorithm
    """
    val = x1+(x2-x1)*(1-y1)/(y2-y1)  
    return val

def interpLin(x1, y1, x2,  y2, xq):
    """
    General linear interpolation/extrapolation between two points.
    """
    return y1 + (xq-x1)*(y2-y1)/(x2-x1)


def find(array, string = None):
    """
    Returns indices and values of non-zero elements. Similar functionality 
    to MATLAB find() function. 

    Copyright (c) 2022 Ahmed Ibrahim
    """
    # making sure that a np array is used
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


def evaluateGrid(xx, yy, charAvg, stdevData, EllipseKFact):
    """
    Evaluates the sampling grid of the marching squares algorithm used to
    extract corridors. At each sampling point, every ellipse in the
    characteristic average is computed, with the maximum value being saved. 

    Parameters
    ----------
    xx: np.ndarray
        Square array of sampling point locations in x
    yy: np.ndarray
        Square array of sampling point locations in y
    charAvg: np.ndarray
        Charactertic average of input signals
    stdevData: np.ndarray
        Standard deviation computed at each point of characteristic average
    EllipseKFact: float
        Factor used to scale the size of the ellipes at each point in the 
        characteristic average, hence scale the extent of the corridors

    Returns
    -------
    zz: np.ndarray
        Square array of sampled values at each sampling point xx, yy

    Notes
    -----
    This function was broken out from main body of arcgen() with the hopes that
    optimization could be performed to improve runtime. Later profiling has 
    shown that this is by no means the slowest part of arcgen(). That would be
    the minimization algorithm and polyxpoly()
    """
    zz = np.zeros(np.shape(xx))   # initalize grid of ellipse values
    nCorr = np.shape(xx)[0]
    # For each grid point, find the max of each standard deviation ellipse
    for iPt in range(nCorr):
        for jPt in range(nCorr):
            zz[iPt,jPt] = np.max((
                                  ((xx[iPt,jPt] - charAvg[:,0])**2 
                                   / (stdevData[:,0]*EllipseKFact)**2 
                                   + (yy[iPt,jPt] - charAvg[:,1])**2 
                                   / (stdevData[:,1]*EllipseKFact)**2)**-1))   
    return zz