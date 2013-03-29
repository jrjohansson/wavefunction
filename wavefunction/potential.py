#
# Potential functions
#
# J Robert Johansson, <robert@riken.jp>
# 

from scipy import *
from nummath import * 



#
# Flux-biased phase qubit potential
#
def U_flux_biased(x, param):

    Ej    = param[0]
    beta  = param[1]
    gamma = param[2]
    
    u = -Ej * (cos(x) - 1 / (2 * beta) * x ** 2 + gamma * x)
    
    return u


#
# Derivative of Current-biased phase qubit potential
#
def dU_dt_flux_biased(x, param):

    if len(x) == 1:
        X = arange(x-0.01,x+0.01,0.001)
        U = U_flux_biased(X, param)
        dU = derivative(U,X)
        du = dU[10]
    else:
        U = U_flux_biased(x, param)
        du = derivative(U,x)     
    
    return du

