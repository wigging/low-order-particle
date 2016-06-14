"""
Function to return the theta (dimensionless temperature) profile at a certain 
dimensionless point (r) in a 1D solid sphere, cylinder, or slab shape.

References: 
1) Recktenwald 2006
2) Bergman, Lavine, Incropera, Dewitt 2011 from Ch. 5, pg.299-304
"""

# Modules and Other Required Functions
#------------------------------------------------------------------------------

import numpy as np
import scipy.special as sp
from funcRoots import roots

# First and Second Terms of the Theta Function
#------------------------------------------------------------------------------
    
def funcCn(root, b):
    """
    First term in the theta function
    root = root from the zeta, Bi equation
    b = shape factor where 2 sphere, 1 cylinder, 0 slab
    """
    if b == 2:
        Cn = 4*(np.sin(root)-root*np.cos(root)) / (2*root-np.sin(2*root))
    elif b == 1:
        Cn = (2/root) * (sp.j1(root) / (sp.j0(root)**2 + sp.j1(root)**2))
    elif b == 0:
        Cn = (4*np.sin(root)) / (2*root + np.sin(2*root))
    return Cn

def funcDn(r, root, b):
    """
    Second term in the theta function
    root = root from the zeta, Bi equation
    b = shape factor where 2 sphere, 1 cylinder, 0 slab
    """
    if b == 2:
        Dn = (1/(root*r)) * np.sin(root*r)
    elif b == 1:
        Dn = sp.j0(root * r)
    elif b == 0:
        Dn = np.cos(root * r)
    return Dn

# Theta Function
#------------------------------------------------------------------------------

def theta(r, b, z, Bi, Fo):
    """
    Dimensionless temperature for analytical solution of 1D transient heat
    conduction for a solid sphere, cylinder, or slab.
    r = dimensionaless length term to evaluate theta, (-)
    b = shape factor where 2 sphere or 1 cylinder or 0 slab, (-)
    z = range of zeta values to evaluate zeta, Bi equation for positive roots
    Bi = Biot number h*L/k, (-)
    Fo = Fourier number alpha*t/L^2, (-)
    """
    
    rts = roots(z, b, Bi)   # positive roots of the zeta, Bi equation
    n = len(rts)            # number of positive roots
    
    # initial dimensionless temperature at first root
    theta = funcCn(rts[0], b)*np.exp(-rts[0]**2 * Fo)*funcDn(r, rts[0], b)
    
    # summation of theta for the remaining roots
    for i in range(1, n):
        dTheta_o = funcCn(rts[i], b)*np.exp(-rts[i]**2 * Fo)*funcDn(r, rts[i], b)
        theta = theta + dTheta_o
        
    return theta    # theta temperature profile evaluated at r
    