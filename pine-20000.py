"""
Compare temperature profiles of 1-D and 3-D models for DF = 20 mm dry loblolly
pine particle. Heat capacity as function of temperature and constant thermal
conductivity. Different equivalent spherical diameters and characteristic lengths
implemented with 1-D model.

Assumptions:
Convection boundary condition at surface.
Symmetry about the center of the solid.
Heat transfer via radiation assumed to be negligable.
Particle does not shrink or expand in size during pyrolysis.

Reference: Wood Handbook 2010
Requirements: Python 3, NumPy, SciPy, Matplotlib, transhc.py
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

As = 2.354e-4   # surface area of Comsol particle, m^2
v = 1.462e-7    # volume of Comsol particle, m^3
H = 0.02027     # height of Comsol particle, m
W = 0.001877    # width of Comsol particle, m
L = 0.005069    # length of Comsol particle, m

# 1D Transient Heat Conduction in Biomass Particle
# -----------------------------------------------------------------------------

# calculate equivalent spherical diameters and characteristic length
ds = (As/np.pi)**(1/2)  # surface area equivalent sphere diameter, m
dv = (6/np.pi*v)**(1/3) # volume equivalent sphere diameter, m
dsv = (dv**3)/(ds**2)   # surface volume equivalent sphere diameter (Sauter), m
dc = v/As               # characteristic length, m
dh = H                  # height as diameter, m
dw = W                  # width as diameter, m
dl = L                  # length as diameter, m

# number of nodes from center of particle (m=0) to surface (m)
m = 1000

# time vector from 0 to max time
tmax = 40.0                     # max time, s
nt = 1000                       # number of time steps
dt = tmax/nt                    # time step, s
t = np.arange(0, tmax+dt, dt)   # time vector, s

# intraparticle temperature array [T] in Kelvin
# row = time step, column = node point from 0 to m
Ts = hc2(ds, x, k, Gb, h, Ti, Tinf, 2, m, t)   # ds case, b = 2 for sphere
Tv = hc2(dv, x, k, Gb, h, Ti, Tinf, 2, m, t)   # dv case, b = 2 for sphere
Tsv = hc2(dsv, x, k, Gb, h, Ti, Tinf, 2, m, t) # dsv case, b = 2 for sphere
Tc = hc2(dc, x, k, Gb, h, Ti, Tinf, 2, m, t)   # dc case, b = 2 for sphere
Th = hc2(dh, x, k, Gb, h, Ti, Tinf, 2, m, t)   # dh case, b = 2 for sphere
Tw = hc2(dw, x, k, Gb, h, Ti, Tinf, 2, m, t)   # dw case, b = 2 for sphere
Tl = hc2(dl, x, k, Gb, h, Ti, Tinf, 2, m, t)   # dl case, b = 2 for sphere

# volume average temperatures
v = vol(ds, m)          # volumes in the sphere
Ts_vol = Tvol(Ts, v)    # ds volume average temperature profile
Tv_vol = Tvol(Tv, v)    # dv volume average temperature profile
Tsv_vol = Tvol(Tsv, v)  # dsv volume average temperature profile
Tc_vol = Tvol(Tc, v)    # dc volume average temperature profile
Th_vol = Tvol(Th, v)    # dh volume average temperature profile
Tw_vol = Tvol(Tw, v)    # dw volume average temperature profile
Tl_vol = Tvol(Tl, v)    # dl volume average temperature profile

# grab data from text file
txtfile = 'comsol/20000tempsPine.txt'
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
py.plot(t, Ts_vol, lw=2, label='Ds')
py.plot(t, Tv_vol, lw=2, label='Dv')
py.plot(t, Tsv_vol, lw=2, label='Dsv')
py.plot(t, Tc_vol, lw=2, label='Dc')
py.plot(t, Th_vol, lw=2, label='Dh')
py.plot(t, Tw_vol, lw=2, label='Dw')
py.plot(t, Tl_vol, lw=2, label='Dl')
py.plot(t2, Tv, 'co', mec='c', lw=2, label='Tv')
py.plot(t2, Tst, 'ms', mec='m', lw=2, label='Tst')
py.plot(t2, Tc, 'kv', mec='k', lw=2, label='Tc')
py.plot(t2, Tl, 'b^', mec='b', lw=2, label='Tl')
py.plot(t2, Tw, 'g<', mec='g', lw=2, label='Tw')
py.plot(t2, Tsa, 'y>', mec='y', lw=2, label='Tsa')
py.axhline(Tinf, c='k', ls='--')
py.ylim(250, 800)
py.xlim(0, tmax)
py.title('1-D (Dsv) vs 3-D Temperatures for DF = 20 mm')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1)
py.grid()
despine()
