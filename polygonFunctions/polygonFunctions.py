import numpy as np
"""
A series of functions pertaining to polygons to be used with ARCGen. These 
functions serve as replacements for the original matlab functions

"""

def signedpolyarea(xy):
    """ 
    Computes the signed area of a polygon

    When positive y is upward, positive area indicates clockwise

    Returns signed area
    """
    xy = np.asarray(xy)
    area = 0
    j = xy.shape[0] - 1

    for i in range(len(xy)):
        area += (xy[j,0] + xy[i,0]) * (xy[j,1] - xy[i,1])
        j = i
    
    area = area/2

    return area

def ispolycw(xy):
    """
    Returns if polygon is clockwise
    """
    return np.sign(signedpolyarea(xy)) > 0

def polyarea(xy):
    """
    Returns the area of a polygon
    """
    return np.abs(signedpolyarea(xy))

