#
# J Robert Johansson, <robert@riken.jp>
#
from scipy import *

#
# Generate matrix representation of Kronecker-delta(i-x,j-y)
#
def kd(nstates, x, y):

    D = zeros((nstates, nstates), dtype=float)
    
    for i in range(0, nstates):
        for j in range(0, nstates):
            
            if (i-x) == (j-y):
                D[i,j] = 1
            else:
                D[i,j] = 0    

    return D

#
# Generate matrix representation of Kronecker-delta product
# kd(a,c)kd(b,d-X)
# where a,b,c,d = -N:N and X is a constant
#
def kdp(nstates, X):

    N = (nstates-1)/2;
    D = zeros((nstates**2, nstates**2), dtype=float)
    k = arange(-N,N+1)
        
    for aidx in range(0,nstates):
        for bidx in range(0,nstates):
            for cidx in range(0,nstates):
                for didx in range(0,nstates):
                    I = (aidx) * nstates + bidx;
                    J = (cidx) * nstates + didx;
                    
                    if (k[aidx] == k[cidx]) & (k[bidx] == k[didx]+X):
                        D[I,J] = 1
                    else:
                        D[I,J] = 0
    
    return D
