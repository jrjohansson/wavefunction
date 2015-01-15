#
# Construct a matrix which represent the schrodinger equation on matrix form
#
# J Robert Johansson, <robert@riken.jp>
#
import numpy as np


#
# kronecker delta, optionally modify so that it also take the boundary
# conditions into account?
#
def mod_kron(N, n, m):
    return (n == m)


#
# Create matrix representing the SE
#
def schrodinger_matrix(xmin, xmax, N, S_param, U_func, U_param):
    """
    legacy api: use assemble_K and assemble_V instead
    """
    h_ = S_param[0]
    mass = S_param[1]

    delta = (xmin-xmax) / N

    M1 = np.zeros((N, N)).astype(float)
    M2 = np.zeros((N, N)).astype(float)

    for n in range(0, N):
        for m in range(0, N):
            x_n = xmin + (xmax-xmin) * n / N
            M1[n, m] = -h_**2/(2*mass*delta**2) * (
                mod_kron(N, n+1, m) - 2 * mod_kron(N, n, m) +
                mod_kron(N, n-1, m))
            M2[n, m] = U_func(x_n, U_param) * mod_kron(N, n, m)

    M = M1 + M2

    return M


def assemble_K(N, k, x_min, x_max, sparse=False):
    """
    Assemble the matrix representation of the kinetic energy contribution
    to the Hamiltonian.

    k = -hbar**2 / 2 m
    """
    dx = (x_min - x_max) / N

    K = np.zeros((N, N)).astype(np.complex)

    for m in range(0, N):
        for n in range(0, N):
            K[m, n] = k / (dx ** 2) * (
                mod_kron(N, m + 1, n) - 2 * mod_kron(N, m, n) +
                mod_kron(N, m - 1, n))

    return K


def assemble_V(N, u, x_min, x_max, sparse=False):
    """
    Assemble the matrix representation of the potential energy contribution
    to the Hamiltonian.
    """
    V = np.zeros((N, N)).astype(np.complex)

    for m in range(N):
        for n in range(N):
            V[m, n] = u[m] * mod_kron(N, m, n)

    return V


def basis_step_evalute(N, u, x):
    """

    """
    return u


def assemble_u_potential(N, u_func, x, args, sparse=False):
    """

    """
    return u_func(x, args)


def wavefunction_norm(x, psi):
    """
    Calculate the norm of the given wavefunction.
    """

    dx = x[1] - x[0]

    return (psi.conj() * psi).sum() * dx


def wavefunction_normalize(x, psi):
    """
    Normalize the given wavefunction.
    """

    return psi / np.sqrt(wavefunction_norm(x, psi))


def expectation_value(x, operator, psi):
    """
    Evaluate the expectation value of a given operator and wavefunction.
    """

    dx = x[1] - x[0]

    return (psi.conj() * operator * psi).sum() * dx


def inner_product(x, psi1, psi2):
    """
    Evaluate the inner product of two wavefunctions, psi1 and psi2, on a space
    described by X1 and X2.
    """

    dx = x[1] - x[0]

    return (psi1.conj() * psi2).sum() * dx


def derivative(x, psi):
    """
    Evaluate the expectation value of a given operator and wavefunction.
    """

    dx = x[1] - x[0]

    N = len(psi)

    dpsi = np.zeros(N, dtype=np.complex)

    def _idx_wrap(M, m):
        return m if m < M else m - M

    for n in range(N):
        dpsi[n] = (psi[_idx_wrap(N, n+1)] - psi[n-1]) / (2 * dx)

    return dpsi
