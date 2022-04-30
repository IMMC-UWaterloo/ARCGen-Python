"""
A series of functions pertaining to polygons to be used with ARCGen. These 
functions serve as replacements for the original matlab functions

"""

import numpy as np
import numpy.typing as npt

def signedpolyarea(xy: npt.ArrayLike):
    """ 
    Computes the signed area of a polygon. When positive y is upward, positive 
    area indicates clockwise. Only suitable for non-complex polygon (no 
    self-intersections)

    Parameters
    ----------
    xy: np.ndarray
        Two column array containing x-y points of the polygon. Does not need
            to be closed
    
    Returns
    -------
    area: float
        Signed area of provided polygon

    Copyright (c) 2022 Devon C. Hartlen
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
    Determines if the inputted polygon is defined clockwise 

    Parameters
    ----------
    xy: np.ndarray
        Two column array containing x-y points of the polygon. Does not need
            to be closed

    Returns
    -------
    ispolycw: bool
        Returns true if polygon is defined clockwise

    Copyright (c) 2022 Devon C. Hartlen
    """
    return np.sign(signedpolyarea(xy)) > 0


def polyarea(xy):
    """
    Determines the area of the inputted polygon. Polygon does not need to be 
    closed. Only suitable for non-complex polygon (no self-intersections)

    Parameters
    ----------
    xy: np.ndarray
        Two column array containing x-y points of the polygon. Does not need
            to be closed
    
    Returns
    -------
    area: float
        Unsigned area of provided polygon

    Copyright (c) 2022 Devon C. Hartlen
    """
    return np.abs(signedpolyarea(xy))


def segintersect(xy1, xy2):
    """
    Determines if line segments xy1, xy2 intersect
    xy1, xy2 are 2x2 np arrays where rows are points, columns are xy

    Determines if two line segments intersect

    Parameters
    ----------
    xy1: np.ndarray
        2x2 array with x coordinate in the first column and y coordinate in 
            the second. Two points are defined by row.
    xy2: np.ndarray
        2x2 array with x coordinate in the first column and y coordinate in 
            the second. Two points are defined by row.

    Returns
    -------
    bool
        Returns True if segments intersect. 
    
    Copyright (c) 2022 Devon C. Hartlen
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


def inpolygon(poly: npt.ArrayLike, point: npt.ArrayLike):
    """
    Determines if the point is in the polygon. Uses ray-casting method. Not
    suitable for self-intersecting polygons, but handles non-convex polygons
    fine. Does not assume closed polygons. 

    Parameters
    ----------
    poly: np.ndarray
        Two column array of points defining a polygon. Polygon must be closed
            (i.e. first and last entry must be the same)
    point: np.ndarray
        1x2 array defining the point to be tested. First column is x coordinate
            second column is y corrdinate. 

    Returns
    -------
    bool
        True if point is inside the polygon
    
    Copyright (c) 2022 Devon C. Hartlen
    """
    poly = np.asarray(poly)
    point = np.asarray(point)
    maxHorizontal = np.max(poly[:,0])

    if poly.shape[0] < 3:
        print('Polygons must have at least 3 points')
        return False

    # Define horizontal ray based on point
    ray = np.vstack((point, [1.1*maxHorizontal, point[1]]))

    # cycle through all line segments,    
    prev = poly.shape[0] - 1
    intcount = 0
    for curr in range(poly.shape[0]):
        # Define current segment to check and test against ray
        trialseg = np.vstack((poly[prev,:], poly[curr,:]))
        if segintersect(trialseg, ray):
            # Check edge case of point collinear with current segment.
            if (_segorientation(poly[prev,:], poly[curr,:], point) == 0):
                # If collinear, assume inside and short circuit counting
                return _segcollinear(poly[prev,:], point, poly[curr,:])
            
            intcount += 1
        prev = curr
    # If the intersect count is odd, then point is in polygon. 
    return (intcount % 2 == 1)


def polyxpoly(poly1: npt.ArrayLike, poly2: npt.ArrayLike):
    """
    Computes the intersection of two polynomials. Does not assumed closed 
    polygons. 

    Parameters
    ----------
    poly1: np.ndarray
        Two column array defining closed polygon (first and last index should
            be the same)
    poly2: np.ndarray
        Two column array defining closed polygon (first and last index should
            be the same)

    Returns
    -------
    interVals: np.ndarray
        Two column array (x, y) of points where the two polygons intersect
    interInds: np.ndarray
        Two column array specifying the indices were the two polygons
            intersect. Column 1 are indices of the first polygon. Returns
            interpolated indices as floats. Use floor() or ceil() to get 
            integer indices on either side of intersection. 

    Copyright (c) 2022 Devon C. Hartlen
    """
    # Loop through each polygon and find intersections
    npoly1 = poly1.shape[0]
    npoly2 = poly2.shape[0]
    indApproxInd = np.zeros_like(poly1)
    count = 0
    for i in range(npoly1-1):
        # Build segment to be tested
        iseg = np.vstack((poly1[i,:], poly1[i+1,:]))
        for j in range(npoly2-1):
            jseg = np.vstack((poly2[j,:], poly2[j+1,:]))
            # Check if segments intersect
            if segintersect(iseg, jseg):
                indApproxInd[count,:] = [i,j]
                count += 1
    indApproxInd = indApproxInd[:count,:].astype(int) # cull
    
    # now that we have intersecting segments, find the actual intecepts
    interVals = np.zeros((count,2))
    interInds = np.zeros((count,2))
    for i,ind in enumerate(indApproxInd):
        # Setu A matrix and B vector for linalg
        A = np.zeros((4,4))
        A[0:2,2] = -1
        A[2:4,3] = -1
        A[0,0] = poly1[ind[0]+1,0] - poly1[ind[0],0]
        A[2,0] = poly1[ind[0]+1,1] - poly1[ind[0],1]
        A[1,1] = poly2[ind[1]+1,0] - poly2[ind[1],0]
        A[3,1] = poly2[ind[1]+1,1] - poly2[ind[1],1]

        B = np.array([
            [-poly1[ind[0],0]],
            [-poly2[ind[1],0]],
            [-poly1[ind[0],1]],
            [-poly2[ind[1],1]],
        ])
        # Solve for both reduced indices and actual intercept locations
        X = np.linalg.solve(A,B)
        # break into indices and locations
        interVals[i,:] = [X[2], X[3]]
        interInds[i,:] = [ind[0]+X[0], ind[1]+X[1]]

    return interVals, interInds




def _segorientation(p1, p2, p3):
    """
    Returns to orientation of three provided points used to compute if segments
    intercept

    Parameters
    ----------
    p1: array-like
        (x,y) location of point 1
    p2: array-like
        (x,y) location of point 2
    p3: array-like
        (x,y) location of point 3

    Returns
    -------
    int
        Either -1, 0, 1 depending on orientation case. 

    See Also
    --------
    segintersect()
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
    Checks if p2 is collinear with p1 and p3.

    Parameters
    ----------
    p1: array-like
        (x,y) location of point 1
    p2: array-like
        (x,y) location of point 2
    p3: array-like
        (x,y) location of point 3

    Returns
    -------
    bool
        True if p2 is collienar with p1 and p3

    See Also
    --------
    segintersect()
    """
    if ( (p2[0] <= max(p1[0], p3[0])) and (p2[0] >= min(p1[0], p3[0]))
        and (p2[1] <= max(p1[1], p3[1])) and (p2[1] >= min(p1[1], p3[1])) ):
        return True
    else:
        return False