"""
Compare volume average temperature profiles, Tv, from 1-D model and 3-D Comsol
simulation of white oak particles with Feret diameters of DF = 200 um to 20 mm.
Surface area to volume diameter, Dsv, is used for the 1-D model.

Requirements: Python 3, NumPy, SciPy, Matplotlib, funcHeatCond, funcOther
"""

import numpy as np
import matplotlib.pyplot as py
from funcHeatCond import hc2
from funcOther import vol, Tvol, dsv

# Parameters
# -----------------------------------------------------------------------------

Gb = 0.72       # basic specific gravity, Wood Handbook Table 4-7, (-)
k = 0.16        # thermal conductivity, W/mK
x = 0           # moisture content, %
h = 350         # heat transfer coefficient, W/m^2*K
Ti = 293        # initial particle temp, K
Tinf = 773      # ambient temp, K

# Comsol Data for Particle Geometry and Temperatures
# -----------------------------------------------------------------------------

# geometry and temperature data for DF = 200 um
sa200 = 5.355e-8                    # surface area of Comsol particle, m^2
v200 = 8.895e-13                    # volume of Comsol particle, m^3
file200 = 'comsol/200tempsOak.txt'  # time and temperatures
t200, Tv200, _, _, _, _, _ = np.loadtxt(file200, skiprows=5, unpack=True)

# geometry and temperature data for DF = 400 um
sa400 = 1.879e-7                    # surface area of Comsol particle, m^2
v400 = 5.553e-12                    # volume of Comsol particle, m^3
file400 = 'comsol/400tempsOak.txt'  # time and temperatures
t400, Tv400, _, _, _, _, _ = np.loadtxt(file400, skiprows=5, unpack=True)

# geometry and temperature data for DF = 700 um
sa700 = 4.836e-7                    # surface area of Comsol particle, m^2
v700 = 2.11e-11                     # volume of Comsol particle, m^3
file700 = 'comsol/700tempsOak.txt'  # time and temperatures
t700, Tv700, _, _, _, _, _ = np.loadtxt(file700, skiprows=5, unpack=True)

# geometry and temperature data for DF = 1400 um
sa1400 = 1.394e-6                       # surface area of Comsol particle, m^2
v1400 = 8.442e-11                       # volume of Comsol particle, m^3
file1400 = 'comsol/1400tempsOak.txt'    # time and temperatures
t1400, Tv1400, _, _, _, _, _ = np.loadtxt(file1400, skiprows=5, unpack=True)

# geometry and temperature data for DF = 2800 um
sa2800 = 4.614e-6                       # surface area of Comsol particle, m^2
v2800 = 4.011e-10                       # volume of Comsol particle, m^3
file2800 = 'comsol/2800tempsOak.txt'    # time and temperatures
t2800, Tv2800, _, _, _, _, _ = np.loadtxt(file2800, skiprows=5, unpack=True)

# geometry and temperature data for DF = 5400 um
sa5400 = 1.716e-5                       # surface area of Comsol particle, m^2
v5400 = 2.877e-9                        # volume of Comsol particle, m^3
file5400 = 'comsol/5400tempsOak.txt'    # time and temperatures
t5400, Tv5400, _, _, _, _, _ = np.loadtxt(file5400, skiprows=5, unpack=True)

# geometry and temperature data for DF = 10000 um
sa10000 = 5.885e-5                      # surface area of Comsol particle, m^2
v10000 = 1.827e-8                       # volume of Comsol particle, m^3
file10000 = 'comsol/10000tempsOak.txt'  # time and temperatures
t10000, Tv10000, _, _, _, _, _ = np.loadtxt(file10000, skiprows=5, unpack=True)

# geometry and temperature data for DF = 20000 um
sa20000 = 2.354e-4                      # surface area of Comsol particle, m^2
v20000 = 1.462e-7                       # volume of Comsol particle, m^3
file20000 = 'comsol/20000tempsOak.txt'  # time and temperatures
t20000, Tv20000, _, _, _, _, _ = np.loadtxt(file20000, skiprows=5, unpack=True)

# 1-D Transient Heat Conduction using Dsv
# -----------------------------------------------------------------------------

# number of nodes from center of particle (m=0) to surface (m)
m = 1000

# time vector from 0 to max time
tmax = 2.0                      # max time, s
nt = 1000                       # number of time steps
dt = tmax/nt                    # time step, s
t = np.arange(0, tmax+dt, dt)   # time vector, s

tmax2 = 20.0                    # max time for large particles, s
t2 = np.arange(0, tmax2+dt, dt) # time vector for large particles, s

# 1-D Transient Heat Conduction for DF = 200 um
# -----------------------------------------------------------------------------

# surface area to volume equivalent sphere diameter Dsv, m
dsv200 = dsv(sa200, v200)

# intraparticle temperature array [T] in Kelvin for Dsv case, b = 2 for sphere
# row = time step, column = node point from 0 to m
Tsv200 = hc2(dsv200, x, k, Gb, h, Ti, Tinf, 2, m, t)

# volume average temperature at each time step
vol200 = vol(dsv200, m)         # volumes in the sphere
Tvol200 = Tvol(Tsv200, vol200)  # Dsv volume average temperature profile

# 1-D Transient Heat Conduction for DF = 400 um
# -----------------------------------------------------------------------------

# surface area to volume equivalent sphere diameter Dsv, m
dsv400 = dsv(sa400, v400)

# intraparticle temperature array [T] in Kelvin for Dsv case, b = 2 for sphere
# row = time step, column = node point from 0 to m
Tsv400 = hc2(dsv400, x, k, Gb, h, Ti, Tinf, 2, m, t)

# volume average temperature at each time step
vol400 = vol(dsv400, m)         # volumes in the sphere
Tvol400 = Tvol(Tsv400, vol400)  # Dsv volume average temperature profile

# 1-D Transient Heat Conduction for DF = 700 um
# -----------------------------------------------------------------------------

# surface area to volume equivalent sphere diameter Dsv, m
dsv700 = dsv(sa700, v700)

# intraparticle temperature array [T] in Kelvin for Dsv case, b = 2 for sphere
# row = time step, column = node point from 0 to m
Tsv700 = hc2(dsv700, x, k, Gb, h, Ti, Tinf, 2, m, t)

# volume average temperature at each time step
vol700 = vol(dsv700, m)          # volumes in the sphere
Tvol700 = Tvol(Tsv700, vol700)   # Dsv volume average temperature profile

# 1-D Transient Heat Conduction for DF = 1400 um
# -----------------------------------------------------------------------------

# surface area to volume equivalent sphere diameter Dsv, m
dsv1400 = dsv(sa1400, v1400)

# intraparticle temperature array [T] in Kelvin for Dsv case, b = 2 for sphere
# row = time step, column = node point from 0 to m
Tsv1400 = hc2(dsv1400, x, k, Gb, h, Ti, Tinf, 2, m, t)

# volume average temperature at each time step
vol1400 = vol(dsv1400, m)           # volumes in the sphere
Tvol1400 = Tvol(Tsv1400, vol1400)   # Dsv volume average temperature profile

# 1-D Transient Heat Conduction for DF = 2800 um
# -----------------------------------------------------------------------------

# surface area to volume equivalent sphere diameter Dsv, m
dsv2800 = dsv(sa2800, v2800)

# intraparticle temperature array [T] in Kelvin for Dsv case, b = 2 for sphere
# row = time step, column = node point from 0 to m
Tsv2800 = hc2(dsv2800, x, k, Gb, h, Ti, Tinf, 2, m, t)

# volume average temperature at each time step
vol2800 = vol(dsv2800, m)           # volumes in the sphere
Tvol2800 = Tvol(Tsv2800, vol2800)   # Dsv volume average temperature profile

# 1-D Transient Heat Conduction for DF = 5400 um
# -----------------------------------------------------------------------------

# surface area to volume equivalent sphere diameter Dsv, m
dsv5400 = dsv(sa5400, v5400)

# intraparticle temperature array [T] in Kelvin for Dsv case, b = 2 for sphere
# row = time step, column = node point from 0 to m
Tsv5400 = hc2(dsv5400, x, k, Gb, h, Ti, Tinf, 2, m, t2)

# volume average temperature at each time step
vol5400 = vol(dsv5400, m)           # volumes in the sphere
Tvol5400 = Tvol(Tsv5400, vol5400)   # Dsv volume average temperature profile

# 1-D Transient Heat Conduction for DF = 10000 um
# -----------------------------------------------------------------------------

# surface area to volume equivalent sphere diameter Dsv, m
dsv10000 = dsv(sa10000, v10000)

# intraparticle temperature array [T] in Kelvin for Dsv case, b = 2 for sphere
# row = time step, column = node point from 0 to m
Tsv10000 = hc2(dsv10000, x, k, Gb, h, Ti, Tinf, 2, m, t2)

# volume average temperature at each time step
vol10000 = vol(dsv10000, m)            # volumes in the sphere
Tvol10000 = Tvol(Tsv10000, vol10000)   # Dsv volume average temperature profile

# 1-D Transient Heat Conduction for DF = 20000 um
# -----------------------------------------------------------------------------

# surface area to volume equivalent sphere diameter Dsv, m
dsv20000 = dsv(sa20000, v20000)

# intraparticle temperature array [T] in Kelvin for Dsv case, b = 2 for sphere
# row = time step, column = node point from 0 to m
Tsv20000 = hc2(dsv20000, x, k, Gb, h, Ti, Tinf, 2, m, t2)

# volume average temperature at each time step
vol20000 = vol(dsv20000, m)            # volumes in the sphere
Tvol20000 = Tvol(Tsv20000, vol20000)   # Dsv volume average temperature profile

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
py.plot(t200, Tv200, 'co', mec='c', mew=2, mfc='none', label='Tv')
py.plot(t400, Tv400, 'co', mec='c', mew=2, mfc='none')
py.plot(t700, Tv700, 'co', mec='c', mew=2, mfc='none')
py.plot(t1400, Tv1400, 'co', mec='c', mew=2, mfc='none')
py.plot(t2800, Tv2800, 'co', mec='c', mew=2, mfc='none')
py.plot(t, Tvol200, 'r', lw=2, label='0.2 mm')
py.plot(t, Tvol400, 'g', lw=2, label ='0.4 mm')
py.plot(t, Tvol700, 'b', lw=2, label='0.7 mm')
py.plot(t, Tvol1400, 'm', lw=2, label='1.4 mm')
py.plot(t, Tvol2800, 'y', lw=2, label='2.8 mm')
py.axhline(Tinf, c='k', ls='--')
py.ylim(250, 800)
py.xlim(0, tmax)
py.title('White Oak with DF = 200-2800 um')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()

py.figure(2)
py.plot(t5400, Tv5400, 'co', mec='c', mew=2, mfc='none', label='Tv')
py.plot(t10000, Tv10000, 'co', mec='c', mew=2, mfc='none')
py.plot(t20000, Tv20000, 'co', mec='c', mew=2, mfc='none')
py.plot(t2, Tvol5400, lw=2, label ='5.4 mm')
py.plot(t2, Tvol10000, lw=2, label='10 mm')
py.plot(t2, Tvol20000, lw=2, label='20 mm')
py.axhline(Tinf, c='k', ls='--')
py.ylim(250, 800)
py.xlim(0, tmax2)
py.title('White Oak with DF = 5.4-20 mm')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()
