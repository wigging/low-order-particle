"""
Function that returns the positive roots of the zeta, Bi equation of the
analytical solution for 1D transient heat conduction of a sphere, cylinder,
or slab shape.

References: 
1) Recktenwald 2006
2) Bergman, Lavine, Incropera, Dewitt 2011 from Ch. 5, pg.299-304
"""

# Modules and Other Required Functions
#------------------------------------------------------------------------------

import numpy as np
import scipy.optimize as op
from funcZeta import funcZetaSph, funcZetaCyl, funcZetaSlab

# Roots Function
#------------------------------------------------------------------------------

def roots(z, b, Bi):
    """
    Returns a list of positive roots from zeta, Bi equation used for the
    analytical solution of 1D transient heat conduction for a solid sphere,
    cylinder, or slab shape.
    z = range of zeta values to test for positive roots
    b = shape factor where 2 sphere, 1 cylinder, 0 slab
    Bi = Biot number h*L/k, (-)
    """
    
    if b==2:
        fz = funcZetaSph(z, Bi)     # evaluate sphere function at z-values
        func = funcZetaSph          # declare sphere function as roots function
    elif b==1:
        fz = funcZetaCyl(z, Bi)     # evaluate cylinder function at z-values
        func = funcZetaCyl          # declare cylinder function as roots function
    elif b==0:
        fz = funcZetaSlab(z, Bi)    # evaluate slab function at the z-values
        func = funcZetaSlab         # declare slab function as roots function
    
    # sign of values in fz as -/+ 1
    sign = np.sign(fz)
    
    # calculate difference between fz values where a non-zero value indicates
    # the location of a possible root of the function
    diff = np.diff(sign)
    
    # b as the shape factor where 2 sphere, 1 cylinder, 0 slab
    # use the np.where function to find array index of value change
    # note that np.where returns a tuple
    if b == 2 or b==0:
        # location of sphere or slab roots
        # find index of only positive value change thus ignoring singularities
        where = np.where(diff>0)
    elif b==1:
        # location of cylinder roots
        where = np.where(diff)
    
    roots = np.zeros(len(where[0])) # setup an empty array to store roots
    
    for i, j in enumerate(where[0]):
        roots[i] = op.brentq(func, z[j], z[j+1], args=(Bi)) # roots for a given interval
        
    return roots    # return a list of the positive roots

