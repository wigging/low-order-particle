"""
Compare temperature profiles of 1-D solid sphere, cylinder, and cube shapes that
are volume equivalent. Note that due to surface area, sphere heats slowest, for
example sphere < cylinder < cube.
"""

import numpy as np
import matplotlib.pyplot as py
from funcHeatCond import hc3
from funcOther import volume, surf, lump, dsv

# Parameters
# -----------------------------------------------------------------------------

Gb = 0.54       # basic specific gravity, Wood Handbook Table 4-7, (-)
cp = 1800       # heat capacity, J/kg*K
k = 0.12        # thermal conductivity, W/mK
x = 0           # moisture content, %
h = 350         # heat transfer coefficient, W/m^2*K
Ti = 293        # initial particle temp, K
Tinf = 773      # ambient temp, K

# Calculations for Volume and Surface Area
# -----------------------------------------------------------------------------

# calculate height and cube side to give same volume as sphere
d = 0.001                   # diameter of sphere, m
ht = 2/3*d                  # cylinder height of equal volume as sphere, m
s = ((np.pi*d**3)/6)**(1/3) # cube side of equal volume as sphere, m

Vsph, Vcyl, Vcube = volume(d, ht, s)    # solid volume of sphere, cylinder, cube
SAsph, SAcyl, SAcube = surf(d, ht, s)   # surface area of sphere, cylinder, cube

xf = 1000   # convert length from m to mm, 1 m = 1000 mm
vf = 1e9    # convert volume from m^3 to mm^3, 1 m^3 = 1e9 mm^3
sf = 1e6    # convert surface area from m^2 to mm^2, 1 m^2 = 1e6 mm^2

print('PARAMETERS in mm '+'-'*60)
print('d_sph = {:.8} \t h_cyl = {:.8} \t a_cube = {:.8}'.format(d*xf, ht*xf, s*xf))

print('VOLUME '+'-'*70)
print('V_sph {:.8} \t V_cyl {:.8} \t V_cube {:.8}'.format(Vsph*vf, Vcyl*vf, Vcube*vf))

print('SURFACE AREA '+'-'*64)
print('SA_sph {:.8} \t SA_cyl {:.8} \t SA_cube {:.8}'.format(SAsph*sf, SAcyl*sf, SAcube*sf))

# 1D Transient Heat Conduction Method
# -----------------------------------------------------------------------------

# number of nodes from center of particle (m=0) to surface (m)
m = 1000

# time vector from 0 to max time
tmax = 4.0                     # max time, s
nt = 400                        # number of time steps
dt = tmax/nt                    # time step, s
t = np.arange(0, tmax+dt, dt)   # time vector, s

# solid sphere, shell, microstructure temperature array [T] in Kelvin
# row = time step, column = node point from 0 (center) to m (surface)
Tsphere = hc3(d, cp, k, Gb, h, Ti, Tinf, 2, m, t)
Tsphere_avg = [np.mean(row) for row in Tsphere]

Tcyl = hc3(d, cp, k, Gb, h, Ti, Tinf, 1, m, t)
Tcyl_avg = [np.mean(row) for row in Tcyl]

Tslab = hc3(s, cp, k, Gb, h, Ti, Tinf, 0, m, t)
Tslab_avg = [np.mean(row) for row in Tslab]

# 1D Transient Heat Conduction Method (Dsv)
# -----------------------------------------------------------------------------

dcyl = dsv(SAcyl, Vcyl)     # cylinder
dcube = dsv(SAcube, Vcube)  # cube

Tsph_dsv = hc3(d, cp, k, Gb, h, Ti, Tinf, 2, m, t)
Tsph_dsv_avg = [np.mean(row) for row in Tsph_dsv]

Tcyl_dsv = hc3(dcyl, cp, k, Gb, h, Ti, Tinf, 2, m, t)
Tcyl_dsv_avg = [np.mean(row) for row in Tcyl_dsv]

Tcube_dsv = hc3(dcube, cp, k, Gb, h, Ti, Tinf, 2, m, t)
Tcube_dsv_avg = [np.mean(row) for row in Tcube_dsv]

# Lumped Capacitance Method
# -----------------------------------------------------------------------------

Tsph_lump = lump(cp, h, k, Gb, t, Vsph, SAsph, Ti, Tinf)    # sphere
Tcyl_lump = lump(cp, h, k, Gb, t, Vcyl, SAcyl, Ti, Tinf)    # cylinder
Tcube_lump = lump(cp, h, k, Gb, t, Vcube, SAcube, Ti, Tinf) # cube

# Plots
# -----------------------------------------------------------------------------

py.ion()
py.close('all')

def despine():
    ax = py.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    py.tick_params(axis='both', bottom='off', top='off', left='off', right='off')

py.figure(1)
py.plot(t, Tsphere_avg, 'o-', markevery=20, mec='b', mew=2, lw=2, label='sphere_1d')
py.plot(t, Tcyl_avg, '^-', markevery=20, mec='g', mew=2, lw=2, label='cylinder_1d')
py.plot(t, Tslab_avg, 's-', markevery=20, mec='r', mew=2, lw=2, label='slab_1d')
py.ylim([250, 800])
py.xlabel('Time (s)')
py.ylabel('Average Temperature (K)')
py.title('1-D method')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()

py.figure(2)
py.plot(t, Tsph_dsv_avg, 'o-', markevery=20, mec='b', mew=2, lw=2, label='sphere_1d')
py.plot(t, Tcyl_dsv_avg, '^-', markevery=20, mec='g', mew=2, lw=2, label='cylinder_1d')
py.plot(t, Tcube_dsv_avg, 's-', markevery=20, mec='r', mew=2, lw=2, label='cube_1d')
py.ylim([250, 800])
py.xlabel('Time (s)')
py.ylabel('Average Temperature (K)')
py.title('1-D method with Dsv')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()

py.figure(3)
py.plot(t, Tsph_lump, 'o-', markevery=20, mec='b', mew=2, lw=2, label='sphere_1d')
py.plot(t, Tcyl_lump, '^-', markevery=20, mec='g', mew=2, lw=2, label='cylinder_1d')
py.plot(t, Tcube_lump, 's-', markevery=20, mec='r', mew=2, lw=2, label='cube_1d')
py.ylim([250, 800])
py.xlabel('Time (s)')
py.ylabel('Average Temperature (K)')
py.title('Lumped method')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()
