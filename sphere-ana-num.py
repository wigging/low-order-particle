"""
Compare 1-D analytical sphere solution to 1-D numerical and 3-D Comsol solutions
for transient heat conduction in solid sphere with constant k and Cp.

Assumptions:
Convection boundary condition at surface.
Symmetry about the center of the solid.
Heat transfer via radiation assumed to be negligable.
Particle does not shrink or expand in size during pyrolysis.

Reference: Wood Handbook 2010
Requirements: Python 3, NumPy, SciPy, Matplotlib, funcHeatCond, funcTheta, funcOther
"""

import numpy as np
import matplotlib.pyplot as py
from funcHeatCond import hc3
from funcTheta import theta
from funcOther import vol, Tvol

# Parameters
# -----------------------------------------------------------------------------

d = 0.001       # diameter of sphere, m
Gb = 0.54       # basic specific gravity, Wood Handbook Table 4-7, (-)
cp = 1800       # heat capacity, J/kg*K
k = 0.12        # thermal conductivity, W/mK
x = 0           # moisture content, %
h = 350         # heat transfer coefficient, W/m^2*K
Ti = 293        # initial particle temp, K
Tinf = 773      # ambient temp, K

# 1D Numerical Solution for Transient Heat Conduction in Solid Sphere
# -----------------------------------------------------------------------------

# number of nodes from center of particle (m=0) to surface (m)
m = 1000

# time vector from 0 to max time
tmax = 4.0                      # max time, s
dt = 0.01                       # time step, s
nt = tmax/dt                    # number of time steps
t = np.arange(0, tmax+dt, dt)   # time vector, s

# intraparticle temperature array [T] in Kelvin
# row = time step, column = node point from 0 (center) to m (surface)
T = hc3(d, cp, k, Gb, h, Ti, Tinf, 2, m, t)
Tavg = [np.mean(row) for row in T]

# volume average temperature at each time step
v = vol(d, m)
Tv = Tvol(T, v)

# 1D Analytical Solution for Transient Heat Conduction in Solid Sphere
# -----------------------------------------------------------------------------

ro = d/2        # radius of sphere (a.k.a outer radius), m
rs = ro/ro      # dimensionless surface radius, (-)
rc = 1e-12/ro   # dimensionless center radius, (-)

z = np.arange(0, 1250, 0.1)         # range to evaluate the zeta, Bi equation
z[0] = 1e-12                        # prevent divide by zero warning

rho = Gb*1000                       # density, kg/m^3
alpha = k/(rho*cp)                  # thermal diffusivity biomass, m^2/s
Bi = (h*ro)/k                       # Biot number, (-)
Fo = (alpha * t) / (ro**2)          # Fourier number, (-)

# surface temperature where ro for outer surface, b=2 for sphere
thetaRo = theta(rs, 2, z, Bi, Fo)   # dimensionless temperature profile
T_o = Tinf + thetaRo*(Ti-Tinf)      # convert theta to temperature in Kelvin, K

# center temperature where r for center, b=2 for sphere
thetaR = theta(rc, 2, z, Bi, Fo)    # dimensionless temperature profile
T_r = Tinf + thetaR*(Ti-Tinf)       # convert theta to temperature in Kelvin, K

# 3D Solid Sphere Temperature Data from Comsol
# -----------------------------------------------------------------------------

sphere = 'comsol/3d-sphere-temps.txt'
t_sphere, Tv_sphere, Tc_sphere, Ts_sphere = np.loadtxt(sphere, skiprows=5, unpack=True)

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
py.plot(t, T_o, 'o', mec='b', mew=2, mfc='none', markevery=5, label='Ts_ana')
py.plot(t, T_r, 's', mec='g', mew=2, mfc='none', markevery=5, label='Tc_ana')
py.plot(t, T[:, -1], 'r-', lw=2, label='Ts_num')
py.plot(t, T[:, 0], 'b-', lw=2, label='Tc_num')
py.axhline(Tinf, c='k', ls='--')
py.ylim(250, 800)
py.xlim(0, tmax)
py.title('1-D Analytical and 1-D Numerical')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()

py.figure(2)
py.plot(t, T_o, 'o', mec='b', mew=2, mfc='none', markevery=5, label='Ts_ana')
py.plot(t, T_r, 's', mec='g', mew=2, mfc='none', markevery=5, label='Tc_ana')
py.plot(t_sphere, Ts_sphere, 'r-', lw=2, label='Ts_3d')
py.plot(t_sphere, Tc_sphere, 'b-', lw=2, label='Tc_3d')
py.axhline(Tinf, c='k', ls='--')
py.ylim(250, 800)
py.xlim(0, tmax)
py.title('1-D Analytical and 3-D Comsol')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()

py.figure(3)
py.plot(t_sphere, Ts_sphere, 'o', mec='r', mew=2, mfc='none', label='Ts_3d' )
py.plot(t_sphere, Tc_sphere, 's', mec='b', mew=2, mfc='none', label='Tc_3d')
py.plot(t_sphere, Tv_sphere, '^', mec='g', mew=2, mfc='none', label='Tv_3d')
py.plot(t, T[:, -1], 'r-', lw=2, label='Ts_1d')
py.plot(t, T[:, 0], 'b-', lw=2, label='Tc_1d')
py.plot(t, Tv, 'g-', lw=2, label='Tv_1d')
#py.plot(t, Tavg, 'y-', lw=2, label='Tavg_1d')
py.axhline(Tinf, c='k', ls='--')
py.ylim(250, 800)
py.xlim(0, tmax)
py.title('1-D Numerical and 3-D Comsol')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()
