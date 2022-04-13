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

def inpolygon(poly, point):
    """
    Determines if the point is in the polygon. Uses ray-casting method. Not
    suitable for self-intersecting polygons, but handles non-convex polygons
    fine
    poly: nx2 np.array of points defining polygon
    point 1x2 np.array defining point of interest
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
    Checks if p2 is collinear to p1 and p3
    p1, p2, p3 are 2 element arrays representing xy points
    """
    if ( (p2[0] <= max(p1[0], p3[0])) and (p2[0] >= min(p1[0], p3[0]))
        and (p2[1] <= max(p1[1], p3[1])) and (p2[1] >= min(p1[1], p3[1])) ):
        return True
    else:
        return False