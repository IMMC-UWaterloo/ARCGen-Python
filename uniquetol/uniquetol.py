from enum import unique
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
    # Algorithm approach: MATLAB mentions a "lexoconographical approach"
    # 1) Sort first column low high
    # 2) loop through first column. find equal entries
    # 3) for equal entries, sort second column
    print('before sort')
    print(a)
    
    # start by sorting array based on first column
    sortInd = np.argsort(a[:,0],axis=0)
    a = a[sortInd,:]

    print('after sort')
    print(a)

    # Now, for each set of matched (exact matched) first column entries, sort 
    # the second column and swap as needed. Do not alter second column
    i = 0
    while i < a.shape[0]:
        index = (a[i:,0] == a[i,0])
        index = np.flatnonzero(index == True)
        if (len(index) > 1):
            subset = a[i+index,:]
            subsetIndex = np.argsort(subset[:,1], axis=0)
            a[i:i+len(index),:] = subset[subsetIndex,:]
            i += len(index)

        else:    
            i += 1
    
    # We now have a lexoconographic sort of a
    print(a)

    # Culling: initialize a "saveindex" that will mark which rows of a are to 
    # be saved. Then, go through a. Look ahead to find any rows within tolerance
    # and mark those to be removed. Only keep first entry. 
    saveIndex = np.ones((a.shape[0]))
    for i in range(a.shape[0]):
        if saveIndex[i] == 0:
            continue
        else:
            # I can probably replace this with vectorized math
            for j in range(i+1, a.shape[0]):
                if (np.max(np.abs( a[i,:] - a[j,:] )) < tol):
                    saveIndex[j] = 0
    
    # Now cull array
    a = a[saveIndex.astype(bool),:]
    print('Culled Array')
    print(a)

    # why do I need this second set of loops? Can I just do my culling in the 
    # first sort?


if __name__ == '__main__':
    a = np.array([
        [1, 1+1e-5],
        [2-1e-5, 2+1e-5],
        [1,1], 
        [2+1e-5, 2],
        [1+1e-5, 2],
        [2, 1],
        [0,0], 
        [1,2+1e-5],
        [2-1e-5, 1],
        [2, 2],
    ])

    # uniquetol(a, 1e-4)
    uniquetol(a, 1e-4)



