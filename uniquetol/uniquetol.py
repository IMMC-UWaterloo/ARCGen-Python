import numpy as np
import numpy.typing as npt

def uniquetol(a:npt.ArrayLike, tol=1e-6):
    """
    UNIQUETOL: returns the unique elements of input array a within tolerance 
    tol. A sorted array is returned. Only nx2 arraylike arrays are suported
    
    Tolerance check is performed between two indices, u & v, per:
    norm(a(u)-a(v)) < tol

    a = nx2 numpy array
    tol = absolute tolerance
    """
    col1 = a[:,0]
    col2 = a[:,1]

    col1sortind = np.argsort(col1)
    col2sortind = np.argsort(col2)

    print(col1sortind)
    print(col2sortind)

    sortedcol1 = col1[col1sortind]
    sortedcol2 = col2[col2sortind]

    difcol1 = np.diff(sortedcol1)
    difcol2 = np.diff(sortedcol2)

    samecol1 = difcol1 < tol
    samecol2 = difcol2 < tol
    
    samecol1 = np.hstack((np.asarray(0), samecol1))
    samecol2 = np.hstack((0, samecol2))
    

    print(samecol1)
    print(samecol2)

    dupIndex = np.ones((samecol1.shape[0]))
    dupIndex = dupIndex.astype(bool)

    for i, index1 in enumerate(col1sortind):
        index2 = np.argwhere(col2sortind == index1)
        if (samecol1[index1] == 1) and (samecol2[index2] == 1):
            dupIndex[i] = 0
            print('row {0} (indices {1}, {2}) is a duplicate'.format(index1, index1, index2))

    print(dupIndex)
    outputarray = a[col1sortind,:]
    print(a)
    outputarray = outputarray[dupIndex,:]

    print(outputarray)



if __name__ == '__main__':
    a = np.array([
        [0,0],
        [1, 1.0000001],
        [1,1], 
        [1.0000001, 2],
        [0.999999, 2],
        [1,2]
    ])

    # uniquetol(a, 1e-4)
    print(np.isclose(a))



