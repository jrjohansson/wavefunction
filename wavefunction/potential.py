#
# Potential functions
#
# J Robert Johansson, <robert@riken.jp>
# 

from scipy import *
from nummath import * 

#
# Harmonic oscillator potential 
#
def U_ho(x, param):

    omega = param[0]
    m     = param[1]
    x0    = param[2]
    
    u = 1/2.0 * m * (omega ** 2) * ((x - x0) ** 2)

    return u

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
# Current-biased phase qubit potential (washboard potential)
#
def U_current_biased(x, param):

	Ej = param[0]
	Ic = param[1] 
	Ib = param[2]

	u = - Ej * ( cos(x) + Ib/Ic * x )

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

