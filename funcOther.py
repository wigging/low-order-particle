"""
Utility functions for project.
"""

import numpy as np

# Functions
# -----------------------------------------------------------------------------

def lump(cp, h, k, Gb, t, V, SA, Ti, Tinf):
    """
    Lumped capacitance method
    """
    Lc = V/SA                   # characteristic length for sphere, m
    Bi = (h*Lc)/k               # Biot number Eq 5.10, (-)

    rho = Gb*1000               # density, kg/m^3
    alpha = k/(rho*cp)          # thermal diffusivity, m^2/s
    Fo = (alpha*t)/(Lc**2)      # Fourier number Eq 5.12, (-)

    phi = np.exp(-Bi*Fo)        # dimensionless temperature Eq. 5.13, (-)
    T = Tinf + (Ti-Tinf)*phi    # temperature, K
    return T


def volume(d, h, s):
    """
    Return volume of a sphere, cylinder, and cube.
    d = diameter of sphere, m
    h = height of cylinder, m
    s = side of cube, m
    """
    sphere = (np.pi*d**3)/6
    cylinder = (np.pi*d**2)/4*h
    cube = s**3
    return sphere, cylinder, cube


def surf(d, h, s):
    """
    Return surface area of a sphere, cylinder, and cube.
    d = diameter of sphere, m
    h = height of cylinder, m
    s = side of cube, m
    """
    sphere = np.pi*d**2
    cylinder = np.pi*d*h + (np.pi*d**2)/2
    cube = 6*s**2
    return sphere, cylinder, cube


def vol(d, m):
    """
    Create a list of volumes that comprise a solid sphere. Volume at center is a
    sphere volume while subsequent volumes are a spherical shell volume. Units
    are determined by diameter parameter d which is assumed in meters for
    documentatation purposes.

    Parameters
    ----------
    d = diameter of sphere, m
    m = number of nodes from center to surface, -

    Returns
    -------
    v = list of volumes, m^3
    """
    r = d/2     # radius of sphere, m
    nr = m-1    # number of radius steps, -
    dr = r/nr   # radius step as delta r, m

    rad = np.linspace(0, r, nr)     # array of each radius step
    vc = (4/3)*np.pi*(dr**3)        # volume of sphere center

    v = []                          # store volumes
    v.append(vc)                    # add sphere center volume

    # calculate each shell volume and store in v
    for i in rad[:-1]:
        ri = i
        ro = i+dr
        vshell = (4/3)*np.pi*(ro**3) - (4/3)*np.pi*(ri**3)
        v.append(vshell)

    # return list of volumes based on number of node points
    return v


def Tvol(T, v):
    """
    Determine volume average temperature for a sphere. First, calculates average
    temperature between each node point to represent the temperature of that
    particular volume. Uses center sphere volume and shell volumes as weights to
    calculate the volume average temperature of entire sphere.

    Parameters
    ----------
    T = array of temperatures from 1-D numerical solution, K
    v = list of volumes inside sphere based on number of node points, m^3

    Returns
    -------
    Tv = volume average temperature of sphere as weighted mean, K
    """
    Tv = [] # store volume average temperatures, K

    # calculate volume average temperature for each row then store in Tv
    for row in T:
        idx = np.arange(0, len(row))    # indices for each temperature in row
        row_avg = []    # store average temperatures between each node point

        for i in idx[:-1]:
            ti = row[i]             # inner temperature
            to = row[i+1]           # outer temperature
            tavg = (ti+to)/2        # average temperature
            row_avg.append(tavg)    # store average temperature

        # store weighted mean of row average temperatures into Tv
        row_wtmean = np.average(row_avg, weights=v)
        Tv.append(row_wtmean)

    # return list of volume average temperature at each time step
    return Tv


def dsv(sa, v):
    """
    Calculate the surface area to volume equivalent spherical diameter as Dsv.
    Units are assumed as meters, m, and square meters, m^2.

    Parameters
    ----------
    sa = surface area, m^2
    v = volume, m^3

    Returns
    -------
    Dsv = surface area to volume diameter, m
    """
    ds = (sa/np.pi)**(1/2)    # surface area equivalent sphere diameter, m
    dv = (6/np.pi*v)**(1/3)   # volume equivalent sphere diameter, m

    # surface area to volume equivalent sphere diameter (Sauter diameter), m
    Dsv = (dv**3)/(ds**2)
    return Dsv


def biot(h, r, k):
    """
    Calculate the dimensionless Biot number represented by Bi.

    Parameters
    ----------
    h = heat transfer coefficient, W/m^2K
    r = radius of particle, m
    k = thermal conductivity, W/mK

    Returns
    -------
    Bi = Biot number, -
    """
    Bi = (h*r) / k
    return Bi


def py1(k, kr, rho, cp, r):
    """
    Calculate the pyrolysis number.

    Parameters
    ----------
    k = thermal conductivity, W/mK
    kr = rate constant, 1/s
    rho = density, kg/m^3
    cp = heat capacity, J/kgK
    r = radius, m

    Returns
    -------
    py = pyrolysis number, -
    """
    py = k / (kr * rho * cp * (r**2))
    return py


def py2(h, kr, rho, cp, r):
    """
    Calculate the pyrolysis number.

    Parameters
    ----------
    h = heat transfer coefficient, W/m^2K
    kr = rate constant, 1/s
    rho = density, kg/m^3
    cp = heat capacity, J/kgK
    r = radius, m

    Returns
    -------
    py = pyrolysis number, -
    """
    py = h / (kr * rho * cp * r)
    return py
