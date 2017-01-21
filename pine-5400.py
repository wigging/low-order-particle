"""
Surface, volume, and center temperatures from 1-D and 3-D models for DF = 5.4 mm
dry loblolly pine particle. Heat capacity as function of temperature and constant
thermal conductivity. Dsv used as equivalent spherical diameter for 1-D model.

Assumptions:
Convection boundary condition at surface.
Symmetry about the center of the solid.
Heat transfer via radiation assumed to be negligable.
Particle does not shrink or expand in size during pyrolysis.

Reference for Cp and k: Wood Handbook 2010
Requirements: Python 3, NumPy, SciPy, Matplotlib, funcHeatCond, funcOther
"""

import numpy as np
import matplotlib.pyplot as py
from funcHeatCond import hc2
from funcOther import vol, Tvol

# Parameters
# -----------------------------------------------------------------------------

Gb = 0.54       # basic specific gravity, Wood Handbook Table 4-7, (-)
k = 0.12        # thermal conductivity, W/mK
x = 0           # moisture content, %
h = 350         # heat transfer coefficient, W/m^2*K
Ti = 293        # initial particle temp, K
Tinf = 773      # ambient temp, K

As = 1.716e-5   # surface area of Comsol particle, m^2
v = 2.877e-9    # volume of Comsol particle, m^3

# 1D Transient Heat Conduction in Biomass Particle
# -----------------------------------------------------------------------------

# calculate equivalent spherical diameters and characteristic length
ds = (As/np.pi)**(1/2)  # surface area equivalent sphere diameter, m
dv = (6/np.pi*v)**(1/3)  # volume equivalent sphere diameter, m
dsv = (dv**3)/(ds**2)   # surface volume equivalent sphere diameter (Sauter), m
dc = v/As               # characteristic length, m

# number of nodes from center of particle (m=0) to surface (m)
m = 1000

# time vector from 0 to max time
tmax = 5.0                      # max time, s
nt = 1000                       # number of time steps
dt = tmax/nt                    # time step, s
t = np.arange(0, tmax+dt, dt)   # time vector, s

# intraparticle temperature array [T] in Kelvin
# row = time step, column = node point from 0 to m
Tsv = hc2(dsv, x, k, Gb, h, Ti, Tinf, 2, m, t) # dsv case, b = 2 for sphere

# volume average temperatures
v = vol(ds, m)          # volumes in the sphere
Tsv_vol = Tvol(Tsv, v)  # dsv volume average temperature profile

# grab data from text file
txtfile = 'comsol/5400tempsPine.txt'
t2, Tv, Tst, Tc, Tl, Tw, Tsa = np.loadtxt(txtfile, skiprows=5, unpack=True)

# Plot Results
# -----------------------------------------------------------------------------

py.ion()
py.close('all')

def despine():
    ax = py.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    py.tick_params(axis='both', bottom='off', top='off', left='off', right='off')

py.figure(1)
py.plot(t, Tsv[:, -1], 'r-', lw=2, label='Ts_1d')
py.plot(t, Tsv_vol, 'g-', lw=2, label='Tv_1d')
py.plot(t, Tsv[:, 0], 'b-', lw=2, label='Tc_1d')
py.plot(t2, Tsa, 'ro', mec='r', mew=2, mfc='none', label='Ts_3d')
py.plot(t2, Tv, 'g^', mec='g', mew=2, mfc='none', lw=2, label='Tv_3d')
py.plot(t2, Tc, 'bs', mec='b', mew=2, mfc='none', label='Tc_3d')
py.axhline(Tinf, c='k', ls='--')
py.ylim(250, 800)
py.xlim(0, tmax)
py.title('Surface, Volume, Center Temperatures for DF = 5.4 mm')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1)
py.grid()
despine()

# plot black and white figure

py.figure(2)
py.plot(t, Tsv[:, -1], c='k', ls='-', lw=2, label='Ts_1d')
py.plot(t, Tsv_vol, c='k', ls='--', lw=2, label='Tv_1d')
py.plot(t, Tsv[:, 0], c='k', ls=':', lw=3, label='Tc_1d')
py.plot(t2, Tsa, ls='', marker='o', markevery=2, mec='k', mew=2, mfc='none', label='Ts_3d')
py.plot(t2, Tv, ls='', marker='^', markevery=2, mec='k', mew=2, mfc='none', label='Tv_3d')
py.plot(t2, Tc, ls='', marker='s', markevery=2, mec='k', mew=2, mfc='none', label='Tc_3d')
py.axhline(Tinf, c='k', ls='--')
py.ylim(250, 800)
py.xlim(0, tmax)
# py.title('Surface, Volume, Center Temperatures for DF = 5.4 mm')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1, frameon=False)
despine()

