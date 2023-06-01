# ARCGen-Python - Arc-length Response Corridor Generation

ARCGen-Python is a port of the original [ARCGen for MATLAB](https://github.com/IMMC-UWaterloo/ARCGen). The python version of ACRGen is not updated as regularly as as MATLAB source code. 

Biofidelity response corridors are commonly used to assess the performance of surrogates such as computational models or anthropomorphic test devices while capturing the variability of experimental data. ARCGen represents a generalized method for computing response corridors and the characteristic average of experimental data capable of accommodating most types of input signals, including experimental data that is time-based, cross-variable, non-monotonic, and/or hysteretic. 

This document provides information on how to use ARCGen as well as a high-level overview of the methodologies that ARCGen uses. Exampels can be found in the 'ExampleCasesAndDatasets' folder. While an overview of ARCGen's operation is provided below, please refer to please refer to [Hartlen and Cronin (2022)](https://www.frontiersin.org/article/10.3389/fbioe.2022.843148) for detailed, rigorous coverage. Examples and datasets provided with this code are taken from Hartlen and Cronin (2022). 

ARCGen-Python is available for Python 3.8 can be installed directly from PyPI as follows. 

```
pip install arcgen-python
```

ARCGen-Python is released under the open-sourced GNU GPL v3 license. No warranty or guarantee of support is provided. The authors hold no responsibility for the validity, accuracy, or applicability of any results obtained from this code.

# Referencing

If you use ARCGen-Python in published research, please use the following citation in your research. 

Hartlen D.C. and Cronin D.S. (2022), "Arc-Length Re-Parametrization and Signal Registration to Determine a Characteristic Average and Statistical Response Corridors of Biomechanical Data." *Frontiers in Bioengineering and Biotechnology* 10:843148. doi: 10.3389/fbioe.2022.843148

Bibtex format:
```
@article{Hartlen_Cronin_2022,
  AUTHOR={Hartlen, Devon C. and Cronin, Duane S.},   
  TITLE={Arc-Length Re-Parametrization and Signal Registration to Determine a Characteristic Average and Statistical Response Corridors of Biomechanical Data},      
  JOURNAL={Frontiers in Bioengineering and Biotechnology},      
  VOLUME={10},      
  YEAR={2022},      
  URL={https://www.frontiersin.org/article/10.3389/fbioe.2022.843148},       
  DOI={10.3389/fbioe.2022.843148},       
  ISSN={2296-4185},   
}
```
# Contributing

If you find discover any bugs or issues, or would like to suggest improvements, please create an [issue on this github repostiory](https://github.com/DCHartlen/ARCGen-Python/issues). You are invited free to submit pull requests to integrate fixes and features directly into ARCGen, although pull requests will be reviewed before integration. 

Anyone is free to fork ARCGen for thier own work, so long as you follow the requirements of the GNU GPL v3 license.

# Usage
ARCGen-Python is called using the function `arcgen()` and returns a the computed characteristic average, inner corridors, outer corridors, as well as processed signal data and debug data. Input signals can be defined in a number of ways. See the discussion of function parameters below. 

A minimum working example of ARCGen is provided below. Additional examples are provided in the ['ExampleCasesAndDatasets' directory](https://github.com/IMMC-UWaterloo/ARCGen-Python/tree/main/ExampleCasesAndDatasets). 

````python
from arcgen import arcgen

charAvg, innCorr, outCorr, processedSignals, debug = arcgen(inputData)
````

## Parameters
`inputSignals`: Contains the signals from which a characteristic average and corridors can be computed. Signals must be two-dimensional (i.e. acceleration-time, force-displacement) but do not need to contain the same number of points or be sampled at the same sampling frequency. There is no limit to the number of input signals ARCGen can accommodate, but runtime will increase as the number of signals increases. `inputSignals` must be defined in one of the following three formats. 
- a string of the path to a CSV file containing all signals stacked columnwise. 
- A list of np.ndarrays, where each entry of the list is a two-column ndarray containing signal points
- A list of dictionaries, where each entry is a dictionary containing the signal points and a signal identifier. Dictionary format is {'data': np.ndarray, 'specId': str}, where 'data' is a two-column array of signal points and str is a string used as the signal identifier. 

The remaining parameters are optional and have defined default values in 'arcgen()' if not used. 

`nResamplePoints`: An integer defining the number of points re-parameterized signals will contain. Increasing this value will provide a smoother characteristic average and corridors. Default: 250. 

`CorridorRes`: An integer defining the number of grid points the marching squares algorithm uses to extract the corridor. This sampling grid is automatically calculated based on the most extreme possible corridor values. Increasing this value will improve how accurately the corridor is extracted. Default: 250. 

`NormalizeSignals`: A character array used to control whether magnitude normalization is performed before arc-length re-parameterization. It is highly recommended that this option is enabled. Input options: 'on' (default), 'off'

`EllipseKFact`: A float used to scale the size of the ellipses defined by the confidence region at each point of the characteristic average. A value of 1.0 creates a corridor that is one standard deviation wide at each point of the characteristic average. However, this only corresponds to approximately 38% of variability for two-dimensional data. The square root of the chi-squared quantile function can be used to define response corridors with respect to p-value. Default: 1.0 (or corridors of plus and minus one standard deviation).

`MinCorridorWidth`: A float value used to enforce a minimum corridor width based on maximum standard deviation. This option can be useful if corridors do not extend to the beginning or end of the characteristic average due to low variability between input signals. However, enforcing a minimum width is not representative of the actual variability of underlying data. Default: 0.0 (does not enforce a minimum corridor width). 

`nWarpCtrlPts`: An integer defining the number of interior control points used for signal registration. Default: 0 (disables signal registration)

`WarpingPenalty`: A float defining the penalty factor used during signal registration. Default: 1e-2. 

`Diagnostics`: A character array used to activate diagnostic plots. Useful for tracing anomalous behaviours. Options: 'off' (default), 'on', 'detailed'.

`resultsToFile`: A boolean flag indicating if `arcgen()` should print results to file. If `True`, `arcgen()` will output the characteristic average and corridors to a csv file in a folder called 'outputs' in the current working directory. `arcgen()` will also output a plot of the inputsignals superimposed with the characteristic average and corridors in the same directory. Default: `False`

## Outputs
`charAvg`: A two-column array containing the computed characteristic average

`innerCorr`: A two-column array defining the inner or lower portion of the corridor

`outerCorr`: A two-column array defining the outer or upper portion of the corridor

`processedSignals`: Contains information about processed signals. Formatted as a list of dictionaries, where each each signal has its own entry in the list and its own dictionary. All dictionaries are defined as {'data': np.ndarray, 'specId': str, 'xMax': float,'xMin': float, 'yMax': float, 'yMin': float, 'maxAlen': float, 'resampled': np.ndarray, 'warpControlPoints': np.ndarray}. 
- data: Four column array contianing original x,y data (cols 0,1), computed arc-length (col 2) and normalized arc-length (col 3)
- specId: signal identifer
- *Max, *Min: maximum and minimum (x,y) values of each signal. 
- maxAlen: total normalized arc-length of the signal
- 'resampled': Three column array of the re-parameterized, registered signal. Col 0: normalized arc-length, Col 1, 2: x,y data
- warpControlPoints': Two column array of control points used during signal registration

`debug`: Dictionary containing information potential useful for debugging. Dictionary is defined as {'charAvg': np.ndarray, 'stDev': np.ndarray, 'preWarpCorrArray': np.ndarray, 'preWarpMeanCorr': float, 'warpedCorrArray': np.ndarray, warpedMeanCorrScore': float} where 
- charAvg: Two column array of the characteristic average (repeated from first positional output)
- stDev: two column array containing the standard devation in x and y for each corresponding point in the characteristic average
- preWarpCorrArray: x and y correlation scores prior to performing signal registration
- preWarpmeanCorr: mean correlation score prior to peforming signal registration
- warpedCorrArray: x and y correlation scores after performing signal registration
- warpedMeanCorrScore: mean correlation score after performing signal registration

# Overview of Operation
For a detailed description of how ARCGen operates, please refer to Hartlen and Cronin (2022). Only a high-level overview is presented here. 

The operation of ARCGen can be broken into three stages: arc-length re-parameterization (from which ARCGen draws its name), signal registration, and statistical analysis. 

## Arc-length Re-parameterization
Arc-length re-parameterization is the critical feature that allows ARCGen to handle input signals that are non-monotonic in both axes (a behaviour called hysteresis). Arc-length provides a convenient means to define input points with respect to a strictly monotonic value inherently tied to the shape of the signal. 

Before computing the arc-length of each curve, all signals are scaled to eliminate issues of differing magnitude between the axes. This scaling is based on the mean extreme values of all input signals, such that the relative magnitude of each scaled signal relative to one another is maintained. These scaled values are only used for arc-length calculations and are not used later in signal registration or statistical analysis. 

Once signals have been normalized, the arc-length of each signal is computed for each point in the signal. The arc-length is then normalized to itself, such that all signals have a normalized arc-length of 1.0. Finally, signals are resampled such that all signals have points at the same normalized arc-lengths. 

## Signal Registration
One of the underlying assumptions of arc-length re-parametrization is that critical features in the signal appear at approximately the same normalized arc-length. However, if said features are not perfectly aligned, an average can skew or smear the resulting characteristic average. Additionally, features such as significant variability or noise can dramatically affect arc-length calculation, changing where critical features occur with respect to normalized arc-length. 

Signal registration (or curve registration) can be applied to help align critical features of signals. This process introduces a warping function for each input signal that subtly alters how each signal is re-sampled with respect to arc-length to align critical features. These warping functions (strictly monotonic, piecewise Hermite cubic splines) are determined for all signals simultaneously by maximizing cross-correlation between all signals. To ensure that warping functions do not produce highly skewed, unrealistic responses, a penalty function is used to limit the amount of warping introduced. 

Signal registration is an optional process and is not needed if input signals are very highly correlated or strictly monotonic. While some experimentation is needed to select the best number of control points, a rule of thumb would be to set the number of interior control points to the number of inflection points expected in the characteristic average. 

## Statistical Analysis
Following arc-length reparameterization, all input signals will have the same number of points at the same normalized arc-length. If signal registration has been performed, the registered points will be used for statistical analysis. Statistical analysis is undertaken in a point-wise fashion at each normalized arc-length. This analysis assumes points are uncorrelated and distributed according to a two-dimensional normal distribution. Based on this assumption, an elliptical confidence region can be constructed at each normalized arc-length. The length of the major and minor axes of these ellipses are proportional to the standard deviation at each arc-length. However, unlike a one-dimensional confidence region, where a plus and minus one standard deviation region will account for 68% of variability, the same two-dimensional elliptical region only accounts for 38%. To control the size of the region based on variance, the quantile function of the chi-squared distribution at the desired variance or p-value can be used to scale the size of the ellipse (optional input option `EllipseKFact`).

The characteristic average of the input signals is defined as the mean value at each normalized arc-length. The response corridors are the envelope of all ellipses. As there is no closed-form way of extracting this envelope, a marching-squares algorithm is used to extract this envelope numerically. Because the envelope is extracted numerically, it is important that the number of resampling points (`nResamplePoints`) are large enough to ensure that ellipses are sufficiently overlapped to provide a smooth, realistic envelope. Similarly, the resolution of the marching squares grid (`CorridorRes`) should be fine enough to capture the shape of the ellipses correctly. This last feature is similar to ensuring that the mesh of a finite element or computational fluid dynamics simulation is fine enough to resolve features. 

# Change Log

## Version 2023.1
Important fix to `polygonFunction.polyxpoly()` as well as fixes to unit tests. 

## Version 2022.1
This is the first production release of ARCGen-Python. There will inevitably be bugs. If you discover any throughout the course of your usage of ARCGen-Python, [please open an issue ticket](https://github.com/IMMC-UWaterloo/ARCGen-Python/issues).
