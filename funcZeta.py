"""
Functions for zeta, Bi equation for sphere, cylinder, and slab cases for use in
the 1D transient heat conduction analytical solution.

sphere:   1-z*cot(z) = Bi         as f(z) = 1-z*cot(z)-Bi
cylinder: z*(J1(z)/J0(z)) = Bi    as f(z) = z*J1(z)-Bi*J0(z)
slab:     z*tan(z) = Bi           as f(z) = z*tan(z)-Bi

Reference:
Bergman, Lavine, Incropera, Dewitt 2011 from Ch. 5, pg.299-304
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import scipy.special as sp

# Functions
#------------------------------------------------------------------------------

def funcZetaSph(z, Bi):
    """
    zeta, Bi function for sphere as f(z) = 1 - z*cot(z) - Bi
    z = zeta values which are later solved for the positive roots for theta
    Bi = Biot number h*L/k, (-)
    """
    f = 1 - z*(1 / np.tan(z)) - Bi
    return f


def funcZetaCyl(z, Bi):
    """
    zeta, Bi function for cylinder as f(z) = z*J1(z) - Bi*J0(z)
    note that J1 and J0 are bessel functions
    z = zeta values which are later solved for the positive roots for theta
    Bi = Biot number h*L/k, (-)
    """
    f = z*sp.j1(z) - Bi*sp.j0(z)
    return f
    

def funcZetaSlab(z, Bi):
    """
    zeta, Bi function for slab as f(z) = z*tan(z) - Bi
    z = zeta values which are later solved for the positive roots for theta
    Bi = Biot number h*L/k, (-)
    """
    f = z*np.tan(z) - Bi
    return f
    
    