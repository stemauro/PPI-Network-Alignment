import numpy as np
from scipy.sparse import csr_matrix

def adjmat(M, idx):
    for ip in idx:
        M[ip[0], ip[1]] = 1
        M[ip[1], ip[0]] = 1

    return csr_matrix(M)

def simmat(M, idx):
    for ip in idx:
        M[ip[0], ip[1]] = 1

    return csr_matrix(M)

def edge_correctness(A, B, X):
    m, n = X.shape
    I = range(m)
    K = range(n)

    preserved = sum(A[i,j] * B[k,l] * X[i,k] * X[j,l] 
                        for l in K for j in I for k in K for i in I)

    return preserved / A.sum()

def main():

    A = np.zeros((6,6), dtype="int32")
    B = np.zeros((8,8), dtype="int32")
    X = np.zeros((6,8), dtype="int32")

    idxa = [(0,1), (2,3), (3,4), (3,5)]
    idxb = [(0,1), (2,3), (3,4), (3,7), (5,6)]
    idxx = [(0,0), (1,1), (2,2)]

    adjA  = adjmat(A, idxa)
    adjB  = adjmat(B, idxb)
    algnX = simmat(X, idxx)
    
    print(adjA)
    print(adjB)
    print(edge_correctness(adjA, adjB, algnX))
    pass



if __name__ == "__main__":
    main()
