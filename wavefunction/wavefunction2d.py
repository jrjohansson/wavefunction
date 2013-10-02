"""

"""

import numpy as np

#
# Fourier series basis for periodic potentials
#

def index_m2v(L, n1, n2):
    return (2*L+1) * (n1+L) + (n2+L)

def index_v2m(L, N):
    
    n1 = np.floor(N / (2 * L + 1)) - L
    n2 = N - (2 * L + 1) * (n1 + L) - L
    
    return n1, n2

def index_s2a(L, s):
    """
    Convert mathematical sum index to an array (vector/matrix) index.
    """
    return s + L

def index_a2s(L, a):
    """
    Convert an array (vector/matrix) index to mathematical sum index.
    """
    
    return a - L 


def delta(a, b):
    return int(a == b)


def convert_m2v(L1, L2, mat):
    
    vec = np.zeros((2*L1 + 1) * (2*L2 + 1), dtype=np.complex)
    
    for m in np.arange(-L1, L1+1):
        for n in np.arange(-L1, L1+1):    
            vec[index_m2v(L1, m, n)] = mat[m + L1, n + L2]
   
    return vec
    
def convert_v2m(L1, L2, vec):

    mat = np.zeros((2 * L1 + 1, 2 * L2 + 1), dtype=np.complex)

    for N in range((2 * L1 + 1) * (2 * L2 + 1)):
        m, n = index_v2m(L1, N)
        mat[index_s2a(L1, m), index_s2a(L2, n)] = vec[N]
    
    return mat



def assemble_K(L1, L2, k11, k12, k22, Tx1, Tx2, sparse=False):
    """
    Assemble the matrix representation of the kinetic energy contribution
    to the Hamiltonian.
    """    
    L1n = 2 * L1 + 1
    L2n = 2 * L2 + 1
    
    K = np.zeros((L1n*L1n, L2n*L2n), dtype=np.complex)
    
    for n1 in np.arange(-L1, L1+1):
        for n2 in np.arange(-L1, L1+1):
            N = index_m2v(L1, n1, n2)
            for m1 in np.arange(-L2, L2+1):
                for m2 in np.arange(-L2, L2+1):
                    M = index_m2v(L2, m1, m2)
                    if delta(n1, m1) * delta(n2, m2):
                        K[N,M] = (k11 * (2 * np.pi * m1 / Tx1) ** 2 + \
                                  k12 * ((2 * np.pi) ** 2 * m1 * m2 / (Tx1 * Tx2)) + \
                                  k22 * (2 * np.pi * m2 / Tx2) ** 2)
                    
    return K


def assemble_V(L1, L2, u, sparse=False):
    """
    Assemble the matrix representation of the potential energy contribution
    to the Hamiltonian.
    """
    L1n = 2 * L1 + 1
    L2n = 2 * L2 + 1
    
    V = np.zeros((L1n*L1n, L2n*L2n), np.complex)
    
    for n1 in np.arange(-L1, L1+1):
        for n2 in np.arange(-L1, L1+1):
            N = index_m2v(L1, n1, n2)
            for m1 in np.arange(-L2, L2+1):
                for m2 in np.arange(-L2, L2+1):
                    M = index_m2v(L2, m1, m2)
                    k1 = n1 - m1
                    k2 = n2 - m2
                    
                    if not (k1 < -L1 or k1 > L1 or k2 < -L2 or k2 > L2):
                        V[N,M] = u[index_s2a(L1, k1), index_s2a(L2, k2)]

    return V


def basis_periodic_evaluate(L1, L2, u, X1, X2, T1, T2):
    """


    """    
    U = np.zeros(X1.shape, dtype=np.complex)

    for n1 in np.arange(-L1, L1+1):
        for n2 in np.arange(-L2, L2+1):            
            U += u[index_s2a(L1, n1), index_s2a(L2, n2)] * np.exp(2j * np.pi * n1 * X1 / T1) * np.exp(2j * np.pi * n2 * X2 / T2)
    
    return U


def assemble_u_periodic_potential(L1, L2, u_func, X1, X2, T1, T2, args, sparse=False):

    globals().update(args)
    
    L1n = 2 * L1 + 1
    L2n = 2 * L2 + 1
    
    U = u_func(X1, X2, args)
    dX1X2 = (X1[0,1] - X1[0,0]) * (X2[1,0] - X2[0,0]) / (T1 * T2)
    
    u = np.zeros((L1n, L2n), dtype=np.complex)

    for n1 in np.arange(-L1, L1+1):
        N1 = n1 + L1
        for n2 in np.arange(-L2, L2+1):
            N2 = n2 + L2
            u[N1, N2] = (U * np.exp(-2j * np.pi * n1 * X1 / T1) * np.exp(-2j * np.pi * n2 * X2 / T2)).sum() * dX1X2
                        
    return u


#
# Step basis for general potentials
#


#
# Position basis functions
#

def wavefunction_norm(X1, X2, psi):
    """
    Calculate the norm of the given wavefunction.
    """    

    dX1, dX2 = X1[0,1] - X1[0,0], X2[1,0] - X2[0,0]
    
    return (psi.conj() * psi).sum() * dX1 * dX2


def wavefunction_normalize(X1, X2, psi):
    """
    Normalize the given wavefunction.
    """    
    
    return psi / np.sqrt(wavefunction_norm(X1, X2, psi))


def expectation_value(X1, X2, operator, psi):
    """
    Evaluate the expectation value of a given operator and wavefunction.
    """    

    dX1, dX2 = X1[0,1] - X1[0,0], X2[1,0] - X2[0,0]
    
    return (psi.conj() * operator * psi).sum() * dX1 * dX2


def inner_product(X1, X2, psi1, psi2):
    """
    Evaluate the inner product of two wavefunctions, psi1 and psi2, on a space
    described by X1 and X2.
    """    

    dX1, dX2 = X1[0,1] - X1[0,0], X2[1,0] - X2[0,0]
    
    return (psi1.conj() * psi2).sum() * dX1 * dX2


def derivative(X1, X2, psi, axis=0):
    """
    Evaluate the expectation value of a given operator and wavefunction.
    """    

    dX1, dX2 = X1[0,1] - X1[0,0], X2[1,0] - X2[0,0]

    M, N = psi.shape
    
    dpsi = np.zeros(psi.shape, dtype=np.complex)

    def _idx_wrap(M, m):
        return m if m < M else m - M

    if axis == 0:
        for m in range(M):
            for n in range(N):
                dpsi[m, n] = (psi[_idx_wrap(M, m+1),n]-psi[m-1,n]) / (2 * dX1)
    else:
        for m in range(M):
            for n in range(N):
                dpsi[m, n] = (psi[m,_idx_wrap(N, n+1)]-psi[m,n-1]) / (2 * dX2)

    return dpsi

