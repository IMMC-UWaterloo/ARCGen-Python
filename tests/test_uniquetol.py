import numpy as np
from arcgen.uniquetol import uniquetol
import pytest

def test_uniquetol_manual():
    """
    Manually defined array, checking two tolerance levels
    """
    testArray = np.array([
        [0, 0],
        [0+1e-6, 0],    # dup, 1 col, 1e-6
        [0-1e-6, 0],    # dup, 1 col, 1e-6
        [0+1e-5, 0],    # dup, 1 col, 1e-5
        [0-1e-5, 0],    # dup, 1 col, 1e-5
        [0, 0+1e-5],    # dup, 2 col, 1e-5
        [0, 0-1e-5],    # dup, 2 col, 1e-5
        [0, 0+1e-6],    # dup, 2 col, 1e-6
        [0, 0-1e-6],    # dup, 2 col, 1e-6
        [0.5, 1],
        [0.5-1e-6, 1],  # Half value edgecase, col 1, 1e-6
        [0.5+1e-6, 1],  # Half value edgecase, col 1, 1e-6
        [0.5-1e-5, 1],  # Half value edgecase, col 1, 1e-5
        [0.5+1e-5, 1],  # Half value edgecase, col 1, 1e-5  
    ])
    # shuffle array before uniquetol
    rng = np.random.default_rng()
    arr = np.arange(testArray.shape[0])
    rng.shuffle(arr)
    testArray = testArray[arr,:]
    uniquerows = uniquetol(testArray, 1e-7)
    print('uniquetol tol=1e-7: {0} entries'.format(uniquerows.shape[0]))
    assert uniquerows.shape[0] == 14

    uniquerows = uniquetol(testArray, 1e-6)
    print('uniquetol tol=1e-6: {0} entries'.format(uniquerows.shape[0]))
    print(uniquerows)
    assert uniquerows.shape[0] == 10

    uniquerows = uniquetol(testArray, 1e-5)
    print('uniquetol tol=1e-5: {0} entries'.format(uniquerows.shape[0]))
    print(uniquerows)
    assert uniquerows.shape[0] == 4

    uniquerows = uniquetol(testArray, 1e-4)
    print('uniquetol tol=1e-4: {0} entries'.format(uniquerows.shape[0]))
    print(uniquerows)
    assert uniquerows.shape[0] == 2

def test_uniquetol_random():
    """
    generate a random 2d array. Add gaussian noise of controlled magnitude. If
    noise limit is below tolerance, all added noise entries will be removed. 
    """
    nhalfarray = 40
    rng = np.random.default_rng(0)
    testArray = rng.random((nhalfarray,2))
    noiseMean = 0
    noiseStDev = 1e-7 # noise is below tolerance of 1e-6
    noiseArray = np.random.default_rng(0).normal(noiseMean, noiseStDev, (nhalfarray,2))

    testArray = 2*testArray-1
    testArray = np.vstack((testArray, testArray+noiseArray))

    uniquerows = uniquetol(testArray, 1e-6)

    print('size before: {0}, size after: {1}'.format(testArray.shape[0], uniquerows.shape[0]))
    assert uniquerows.shape[0] == nhalfarray
    

if __name__ == '__main__':
    test_uniquetol_manual()
    test_uniquetol_random()