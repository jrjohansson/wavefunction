#
# J Robert Johansson, <robert@riken.jp>
#
from scipy import *
from scipy import linalg
from misc import *
from copy import *

#
# Numerical derivative
#
def derivative(f, x):
    """
    Calculate the numerical derivative of f given a discetization of the
    coordinate x
    """
    L = len(x)

    ddx = zeros(L).astype(type(x))

    ddx[0] =  (f[1]-f[0]) / (x[1] - x[0])
 
    for n in range(1,L-1):
        ddx[n] = (f[n+1]-f[n-1]) / (x[n+1] - x[n-1])
   
        ddx[L-1] = (f[L-1]-f[L-2]) / (x[L-1] - x[L-2])
    
        return ddx

#
# Sorted eigenvectors
#
def eigenvectors_sorted(H):

    vv, SS = linalg.eig(H)
        
    if max(abs(imag(vv))) > 1e-6:
       raise Exception('Complex eigenvalues of Hamiltonian.')
    else:
        vv = real(vv)

    evs = map(None, copy(vv), copy(transpose(SS)))
    evs.sort()
    
    for i in range(0,len(vv)):
        vv[i]   = evs[i][0]
        SS[i,:] = evs[i][1]
        
    return vv, transpose(SS)

#
# Sorted eigenvectors
#
def eigenvectors_sorted_new(H):

    evals, evecs = linalg.eig(H)
        
    if max(abs(imag(evals))) > 1e-6:
       raise Exception('Complex eigenvalues of Hamiltonian.')
    else:
       evals = real(evals)

    #evs = map(None, copy(vv), copy(transpose(SS)))
    #evs.sort()

    evs = list(zip(evals, range(len(evals))))
    evs.sort()
    evals, perm = list(zip(*evs))
    
    print "evals =", evals

    print "perm =", perm

    #for i in range(0,len(vv)):
    #    vv[i]   = evs[i][0]
    #    SS[i,:] = evs[i][1]
    #    
    #return vv, transpose(SS)

    evecs_sorted = array([evecs[:,k] for k in perm])

    return evals, evecs_sorted

#
# Sorted eigenvalues
#
def eigenvalues_sorted(H):

    vv, SS = linalg.eig(H)
       
    vv = sort(real(vv))
       
    return vv
