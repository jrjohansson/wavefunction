#
#
#

from scipy import *

#
# Create M-matrix
#
def recastMtoT(nstates, a, b):

    M = zeros((nstates**2,nstates**2), dtype=float)

    for i in range(0,nstates):
       for j in range(0,nstates):
            for p in range(0,nstates):
                for q in range(0,nstates):
                    I = (i) * nstates + j;
                    J = (p) * nstates + q;
                    M[I,J] = a[i,p]*b[q,j];
    
    return M

#
# Create M-matrix
#
def recastMtoV(nstates, M):
    v = zeros(nstates**2, dtype=float)
    for i in range(0,nstates):
        for j in range(0,nstates):
            v[i*nstates+j] = M[i,j]
    return v

#
# Create M-matrix
#
def recastVtoM(nstates, v):

    M = zeros((nstates,nstates), dtype=float)
    for i in range(0,nstates):
        for j in range(0,nstates):
             M[i,j] = v[i*nstates+j] 
        
    return M
