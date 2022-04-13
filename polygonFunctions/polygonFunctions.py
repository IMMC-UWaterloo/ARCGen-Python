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

def segintersect(xy1, xy2):
    """
    Determines if line segments xy1, xy2 intersect
    xy1, xy2 are 2x2 np arrays where rows are points, columns are xy
    """

    # Compute the four required orientations
    orient1 = _segorientation(xy1[0,:], xy1[1,:], xy2[0,:])
    orient2 = _segorientation(xy1[0,:], xy1[1,:], xy2[1,:])
    orient3 = _segorientation(xy2[0,:], xy2[1,:], xy1[0,:])
    orient4 = _segorientation(xy2[0,:], xy2[1,:], xy1[1,:])

    # General case: Intersect if (orient1 != orient2) AND (orient3 != orient4)
    if ((orient1 != orient2) and (orient3 != orient4)):
        # Segments intersect
        return True
    
    # Check special collinear cases
    if ( (orient1 == 0) and _segcollinear(xy1[0,:], xy2[0,:], xy1[1,:])):
        return True

    if ( (orient2 == 0) and _segcollinear(xy1[0,:], xy2[1,:], xy1[1,:])):
        return True

    if ( (orient3 == 0) and _segcollinear(xy2[0,:], xy1[0,:], xy2[1,:])):
        return True

    if ( (orient4 == 0) and _segcollinear(xy2[0,:], xy1[1,:], xy2[1,:])):
        return True

    # otherwise return false
    return False


def _segorientation(p1, p2, p3):
    """
    Returns to orientation of three provided points
    p1, p2, p3 are 2 element arrays representing xy point location
    """

    orient = ( (p2[1]-p1[1]) * ( p3[0]-p2[0]) 
                - (p3[1]-p2[1]) * (p2[0]-p1[0]))
    
    if (orient > 0):
        # clockwise
        return 1
    elif (orient < 0):
        # counter-clockwise
        return -1
    else:
        # collinear case
        return 0

def _segcollinear(p1, p2, p3):
    """
    Checks if p3 is collinear to p1 and p2
    p1, p2, p3 are 2 element arrays representing xy points
    """
    if ( (p2[0] <= max(p1[0], p3[0])) and (p2[0] >= max(p1[0], p3[0]))
        and (p2[1] <= max(p1[1], p3[1])) and (p2[1] >= max(p1[1], p3[1])) ):
        return True
    else:
        return False