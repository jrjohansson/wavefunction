#
# Construct a matrix which represent the schrodinger equation on matrix form
# 
# J Robert Johansson, <robert@riken.jp>
#

from scipy import *

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

    h_   = S_param[0]
    mass = S_param[1]

    delta = (xmin-xmax) / N
  
    M1 = zeros((N,N)).astype(float)
    M2 = zeros((N,N)).astype(float)
  
    for n in range(0, N):
        for m in range(0,N):
            x_n = xmin + (xmax-xmin) * n / N
            M1[n,m] = -h_**2/(2*mass*delta**2) * (mod_kron(N, n+1, m) - 2 * mod_kron(N, n, m) + mod_kron(N, n-1, m));     
            M2[n,m] = U_func(x_n, U_param) * mod_kron(N, n, m);
        
    M = M1 + M2;
    
    return M
