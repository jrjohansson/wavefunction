"""
Utility functions for wavefunction calculations.

J Robert Johansson, <robert@riken.jp>
"""

import numpy as np

def print_matrix(A):
    """
    Print real part of matrix matrix to stdout
    """
    M, N = A.shape

    for m in range(0, M):
        for n in range(0, N):
            val = A[m,n].real
            if val > 0.0:
                print(" %.3f  " % val),
            else:
                print("%.3f  " % val),

        print


def solve_eigenproblem(H):
    """
    Solve an eigenproblem and return the eigenvalues and eigenvectors.
    """
    vals, vecs = np.linalg.eig(H)
    idx = np.real(vals).argsort()
    vals = vals[idx]
    vecs = vecs.T[idx]

    return vals, vecs
