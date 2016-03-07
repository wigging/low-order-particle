"""
Compare temperature profiles of 3-D cube and 3-D sphere in Comsol to 1-D cube and
1-D sphere model. Compare 3-D cube and 3-D sphere Comsol temperature profiles.
Sphere and cube are volume equivalent where sphere diameter is 1 mm and cube
side is 0.806 mm.
"""

import numpy as np
import matplotlib.pyplot as py
from funcHeatCond import hc3
from funcOther import volume, surf, dsv, vol, Tvol

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

# Compare Surface Area and Volume for Sphere, Cylinder, and Cube
# -----------------------------------------------------------------------------

# calculate cylinder height and cube side to give same volume as sphere
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

# 1-D Transient Heat Conduction in Solid Sphere
# -----------------------------------------------------------------------------

# number of nodes from center of particle (m=0) to surface (m)
m = 1000

# time vector from 0 to max time
tmax = 4.0                      # max time, s
nt = 400                        # number of time steps
dt = tmax/nt                    # time step, s
t = np.arange(0, tmax+dt, dt)   # time vector, s

# temperatures using d
# row = time step, column = node point from 0 (center) to m (surface)
T_sphere1d = hc3(d, cp, k, Gb, h, Ti, Tinf, 2, m, t)    # array of temperatues, K

# volume average temperature at each time step based on d
vs = vol(d, m)
Tv_sphere1d = Tvol(T_sphere1d, vs)

# temperatures using dsv
# row = time step, column = node point from 0 (center) to m (surface)
dsv_sphere = dsv(SAsph, Vsph)
T_sphere1dsv = hc3(dsv_sphere, cp, k, Gb, h, Ti, Tinf, 2, m, t)

vs_dsv = vol(dsv_sphere, m)
Tv_sphere1dsv = Tvol(T_sphere1dsv, vs_dsv)

# 1-D Transient Heat Conduction Method in Solid Cube
# -----------------------------------------------------------------------------

# volume average temperatures
dsv_cube = dsv(SAcube, Vcube)  # cube
Tcube_dsv = hc3(dsv_cube, cp, k, Gb, h, Ti, Tinf, 2, m, t)

vc_dsv = vol(dsv_cube, m)
Tv_cube1dsv = Tvol(Tcube_dsv, vc_dsv)

# 3-D Cube and 3-D Sphere Temperature Data from Comsol
# -----------------------------------------------------------------------------

cube = 'comsol/3d-cube-temps.txt'
t_cube3d, Tv_cube3d, Tc_cube3d, Tcnr_cube3d, Ts_cube3d = np.loadtxt(cube, skiprows=5, unpack=True)

sphere = 'comsol/3d-sphere-temps.txt'
t_sphere3d, Tv_sphere3d, Tc_sphere3d, Ts_sphere3d = np.loadtxt(sphere, skiprows=5, unpack=True)

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
py.plot(t_sphere3d, Tv_sphere3d, 'g', lw=2, label='Tv_sphere3d')
py.plot(t, Tv_sphere1d, 'bo', mec='b', mew=2, markevery=20, label='Tv_sphere1d')
py.plot(t, Tv_sphere1dsv, 's', markevery=20, mfc='none', mec='r', mew=2, label='Tv_sphere1dsv')
py.plot(t_cube3d, Tv_cube3d, 'c', lw=2, label='Tv_cube3d')
py.plot(t, Tv_cube1dsv, 's', markevery=20, mfc='none', mec='m', mew=2, label='Tv_cube1dsv')
py.xlim([0, 4])
py.ylim([250, 800])
py.xlabel('Time (s)')
py.ylabel('Average Temperature (K)')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()

py.figure(2)
py.plot(t_cube3d, Tv_cube3d, lw=2, label='Tv')
py.plot(t_cube3d, Tc_cube3d, lw=2, label='Tc')
py.plot(t_cube3d, Tcnr_cube3d, lw=2, label='Tcn')
py.plot(t_cube3d, Ts_cube3d, lw=2, label='Ts')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title('3-D Cube Temperatures, s = 0.806mm')
py.legend(loc='lower right', numpoints=1, frameon=False)
py.grid()
despine()

py.figure(3)
py.plot(t_sphere3d, Tv_sphere3d, lw=2, label='Tv')
py.plot(t_sphere3d, Tc_sphere3d, lw=2, label='Tc')
py.plot(t_sphere3d, Ts_sphere3d, lw=2, label='Ts')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title('3-D Sphere Temperatures, d = 1mm')
py.legend(loc='lower right', numpoints=1, frameon=False)
py.grid()
despine()

py.figure(4)
py.plot(t_cube3d, Ts_cube3d, 'r-', lw=2, label='Ts_cube_3d')
py.plot(t_cube3d, Tc_cube3d, 'r--', lw=2, label='Tc_cube_3d')
py.plot(t_sphere3d, Ts_sphere3d, 'b-', lw=2, label='Ts_sphere_3d')
py.plot(t_sphere3d, Tc_sphere3d, 'b--', lw=2, label='Tc_sphere_3d')
py.xlim([0, 4])
py.ylim([250, 800])
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title('Compare 3-D Cube and Sphere')
py.legend(loc='lower right', numpoints=1, frameon=False)
py.grid()
despine()
