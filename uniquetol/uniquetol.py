
import numpy as np
import numpy.typing as npt
import timeit

def uniquetol(a:npt.ArrayLike, tol=1e-6):
    """
    UNIQUETOL: returns the unique elements of input array a within tolerance 
    tol. A sorted array is returned. Only nx2 arraylike arrays are suported
    
    Tolerance check is performed between two indices, u & v, per:
    norm(a(u)-a(v)) < tol

    a = nx2 numpy array
    tol = absolute tolerance
    """
        # Algorithm approach: MATLAB mentions a "lexoconographical approach"
    # 1) Sort first column low high
    # 2) for each row, find all subsequent rows within tolerance and cull
    # This is an expensive appraoch, but robust. 
    
    # start by sorting array based on first column
    sortInd = np.argsort(a[:,0],axis=0)
    a = a[sortInd,:]

    # Going row-by-row, search ahead for rows within tolerance
    # Probably a way to window this method for faster execution.
    saveindex = np.ones((a.shape[0]))
    for i in range(a.shape[0]):
        # if row is already culled, skip
        if saveindex[i] == 0:
            continue
        # Find all rows within maximum tolerance
        else:
            # within tol is when difference in both columns is below tol
            withintol = np.max(np.abs(a[i,:] - a[i+1:, :]) ,axis=1) < tol
            withintol = np.hstack((np.zeros(i+1).astype(bool), withintol))
            saveindex[withintol.astype(bool)] = 0

    a = a[saveindex.astype(bool)]

    return a


if __name__ == '__main__':
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

    # run uniquetol
    uniquerows = uniquetol(testArray, 1e-5)
    print('uniquetol result: {0} entries'.format(uniquerows.shape[0]))
    print(uniquerows)

    SETUPCODE="""
import numpy as np
from __main__ import uniquetol
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
    """
    TESTCODE="""
uniquetol(testArray, 1e-5)
    """
    v1time = timeit.timeit(
        setup=SETUPCODE,
        stmt=TESTCODE,
        number=10000
    )
    print('uniquetol time for 10000 ={0}'.format(v1time))
