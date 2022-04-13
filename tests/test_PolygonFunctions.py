import polygonFunctions as poly
import pytest

import numpy as np

def test_signedpolyarea():
    """
    Unit test to ensure that signedArea works for normal square
    """
    # Define a clockwise 1x1 square
    xycw = np.array([[0, 0], 
                   [0, 1],
                   [1, 1],
                   [1, 0]])
    area = poly.signedpolyarea(xycw)
    # Test if the sign is positive for cw
    assert np.sign(area) == 1
    # Test that area is 1
    assert np.abs(area) == 1
    print('CW Area = %f' %area)
    
    # Define a counter-clockwise 1x1 square
    xyccw = np.array([[0, 0],
                      [1, 0],
                      [1, 1],
                      [0, 1]])
    area = poly.signedpolyarea(xyccw)
    # Test if sign is negative for ccw
    assert np.sign(area) == -1
    # Ensure that area is 1
    assert np.abs(area) == 1
    print('CCW Area = %f' %area)

def test_ispolycw():
    """
    Test to ensure that polygons are correctly recognized as cw with complex
    polygon
    """
    xyccw = np.array([
        [0,0],
        [2,0],
        [2,1], 
        [1,1],
        [1,2],
        [2,2],
        [2,3],
        [3,3],
        [0,3]
    ])

    xycw = np.flipud(xyccw)

    assert poly.ispolycw(xycw)

    assert poly.ispolycw(xyccw) == 0

def test_polyarea():
    """
    Test area of complex polygon that is CCW defined
    """
    xyccw = np.array([
        [0,0],
        [2,0],
        [2,1], 
        [1,1],
        [1,2],
        [2,2],
        [2,3],
        [3,3],
        [0,3]
    ])

    area = poly.polyarea(xyccw)

    print('Computed Polygon Area %f' %area)

    assert pytest.approx(area) == 5.0

def test_segintersect():
    """
    Test line segment intersections
    """
    # Define simple intersecting lines
    xy1 = np.array([
        [0, 0],
        [1, 1]
    ])
    xy2 = np.array([
        [0, 1],
        [1, 0]
    ])
    intersect = poly.segintersect(xy1, xy2)
    print('case 1: intersecting = {0}'.format(intersect))
    assert (intersect == True)

    # Define simple intersecting lines, but ccw
    xy1 = np.array([
        [0, 0],
        [1, 1]
    ])
    xy2 = np.array([
        [1, 0],
        [0, 1]
    ])
    intersect = poly.segintersect(xy1, xy2)
    print('case 2: intersecting = {0}'.format(intersect))
    assert (intersect == True)

    # define non-intersecting
    xy1 = np.array([
        [0, 0],
        [0, 0]
    ])
    xy2 = np.array([
        [1, 1],
        [1, 1]
    ])
    intersect = poly.segintersect(xy1, xy2)
    print('case 3: collinear, non-intersecting = {0}'.format(intersect))
    assert (intersect == False)

    # define on-segment intersecting
    xy1 = np.array([
        [0, 0],
        [1, 1]
    ])
    xy2 = np.array([
        [1, 0],
        [1, 1]
    ])
    intersect = poly.segintersect(xy1, xy2)
    print('case 4: on-segment intersecting = {0}'.format(intersect))
    assert (intersect == True)


if __name__ == '__main__':
    # test_signedpolyarea()
    # test_ispolycw()
    # test_polyarea()
    test_segintersect()
