"""
Compare 1-D transient heat conduction model to Sadhukhan 2009 Figure 2 cylinder.
"""

import numpy as np
import matplotlib.pyplot as py
from funcHeatCond import hc
from funcOther import dsv
from funcKinetics import kn

# Parameters
# -----------------------------------------------------------------------------

d = 0.02        # wood particle diameter, m
rhow = 682      # density of wood, kg/m^3
Ti = 285        # initial particle temp, K
Tinf = 683      # ambient temp, K
h = 48          # heat transfer coefficient, W/m^2*K
H = -240000     # heat of reaction, J/kg where (-)=exothermic, (+)=endothermic

nt = 2000                       # number of time steps
tmax = 800                      # max time, s
dt = tmax/nt                    # time step, s
t = np.arange(0, tmax+dt, dt)   # time vector

nr = 999    # number or radius steps
r = d/2     # radius of particle, m
dr = r/nr   # radius step, delta r
m = nr+1    # nodes from center m=0 to surface m=steps+1

# 1-D Transient Heat Conduction with Dsv and b = 2 with H
# -----------------------------------------------------------------------------

ht = 0.1                    # height of cylinder, m
v = (np.pi*(d**2)*ht)/4     # volume of cylinder, m^3
sa = np.pi*d*((d/2)+ht)     # surface area of cylinder, m^2
dsv = dsv(sa, v)            # Sauter diameter, m

r_dsv = dsv/2     # radius of particle, m
dr_dsv = r_dsv/nr   # radius step, delta r

TdsvH = np.zeros((len(t), m))
TdsvH[0] = Ti

# density array
# rows = time step, columns = node points from center to surface
pw_dsvH = np.zeros((len(t), m))      # create array for wood density
pc_dsvH = np.zeros((len(t), m))      # create array for char density
pg_dsvH = np.zeros((len(t), m))      # create array for gas density

pw_dsvH[0] = rhow                 # initial wood density at all nodes

# mass fraction array
B_dsvH = np.ones((len(t), m))
C1_dsvH = np.zeros((len(t), m))
C2_dsvH = np.zeros((len(t), m))

# mass fraction vector
# columns = average mass fraction of entire solid at a time step
Ys_dsvH = np.ones(len(t))   # create row vector for mass fraction, Ys=1 for all wood

Yw_dsvH = pw_dsvH[0]/rhow     # wood fraction, Yw=1 all wood, Yw=0 all char

cpw_dsvH = 1112.0 + 4.85 * (TdsvH[0] - 273.15)   # wood heat capacity, J/(kg*K)
kw_dsvH = 0.13 + (3e-4) * (TdsvH[0] - 273.15)    # wood thermal conductivity, W/(m*K)
cpc_dsvH = 1003.2 + 2.09 * (TdsvH[0] - 273.15)   # char heat capacity, J/(kg*K)
kc_dsvH = 0.08 - (1e-4) * (TdsvH[0] - 273.15)    # char thermal conductivity, W/(m*K)

cpbar_dsvH = Yw_dsvH*cpw_dsvH + (1-Yw_dsvH)*cpc_dsvH    # effective heat capacity
kbar_dsvH = Yw_dsvH*kw_dsvH + (1-Yw_dsvH)*kc_dsvH       # effective thermal conductivity
pbar_dsvH = pw_dsvH[0] + pc_dsvH[0]                     # effective density

g_dsvH = np.ones(m)*(1e-10)  # assume initial heat generation is negligible

# solve system of equations [A]{T}={C} where T = A\C for each time step
for i in range(1, nt+1):

    # heat conduction
    TdsvH[i] = hc(m, dr_dsv, 2, dt, h, Tinf, g_dsvH, TdsvH, i, r_dsv, pbar_dsvH, cpbar_dsvH, kbar_dsvH)

    # kinetic reactions
    B_dsvH[i], C1_dsvH[i], C2_dsvH[i], g_dsvH = kn(TdsvH, B_dsvH, C1_dsvH, C2_dsvH, rhow, dt, i, H)

    # update thermal properties
    cpw_dsvH = 1112.0 + 4.85 * (TdsvH[i] - 273.15)
    kw_dsvH = 0.13 + (3e-4) * (TdsvH[i] - 273.15)
    cpc_dsvH = 1003.2 + 2.09 * (TdsvH[i] - 273.15)
    kc_dsvH = 0.08 - (1e-4) * (TdsvH[i] - 273.15)

    # update wood and char density
    pw_dsvH[i] = B_dsvH[i]*rhow
    pc_dsvH[i] = (C1_dsvH[i]+C2_dsvH[i])*rhow

    # update mass fraction vector
    Yw_dsvH = pw_dsvH[i] / (pw_dsvH[i] + pc_dsvH[i])
    cpbar_dsvH = Yw_dsvH*cpw_dsvH + (1-Yw_dsvH)*cpc_dsvH
    kbar_dsvH = Yw_dsvH*kw_dsvH + (1-Yw_dsvH)*kc_dsvH
    pbar_dsvH = pw_dsvH[i] + pc_dsvH[i]
    Ys_dsvH[i] = np.mean(B_dsvH[i] + C1_dsvH[i] + C2_dsvH[i])

# 1-D Transient Heat Conduction with Dsv and b = 2 with H = 0
# -----------------------------------------------------------------------------

Tdsv = np.zeros((len(t), m))
Tdsv[0] = Ti

# density array
# rows = time step, columns = node points from center to surface
pw_dsv = np.zeros((len(t), m))      # create array for wood density
pc_dsv = np.zeros((len(t), m))      # create array for char density
pg_dsv = np.zeros((len(t), m))      # create array for gas density

pw_dsv[0] = rhow                 # initial wood density at all nodes

# mass fraction array
B_dsv = np.ones((len(t), m))
C1_dsv = np.zeros((len(t), m))
C2_dsv = np.zeros((len(t), m))

# mass fraction vector
# columns = average mass fraction of entire solid at a time step
Ys_dsv = np.ones(len(t))   # create row vector for mass fraction, Ys=1 for all wood

Yw_dsv = pw_dsv[0]/rhow     # wood fraction, Yw=1 all wood, Yw=0 all char

cpw_dsv = 1112.0 + 4.85 * (Tdsv[0] - 273.15)   # wood heat capacity, J/(kg*K)
kw_dsv = 0.13 + (3e-4) * (Tdsv[0] - 273.15)    # wood thermal conductivity, W/(m*K)
cpc_dsv = 1003.2 + 2.09 * (Tdsv[0] - 273.15)   # char heat capacity, J/(kg*K)
kc_dsv = 0.08 - (1e-4) * (Tdsv[0] - 273.15)    # char thermal conductivity, W/(m*K)

cpbar_dsv = Yw_dsv*cpw_dsv + (1-Yw_dsv)*cpc_dsv     # effective heat capacity
kbar_dsv = Yw_dsv*kw_dsv + (1-Yw_dsv)*kc_dsv        # effective thermal conductivity
pbar_dsv = pw_dsv[0] + pc_dsv[0]                    # effective density

g_dsv = np.ones(m)*(1e-10)  # assume initial heat generation is negligible

# solve system of equations [A]{T}={C} where T = A\C for each time step
for i in range(1, nt+1):

    # heat conduction
    Tdsv[i] = hc(m, dr_dsv, 2, dt, h, Tinf, g_dsv , Tdsv, i, r_dsv, pbar_dsv, cpbar_dsv, kbar_dsv)

    # kinetic reactions
    B_dsv[i], C1_dsv[i], C2_dsv[i], g_dsv = kn(Tdsv, B_dsv, C1_dsv, C2_dsv, rhow, dt, i, 0)

    # update thermal properties
    cpw_dsv = 1112.0 + 4.85 * (Tdsv[i] - 273.15)
    kw_dsv = 0.13 + (3e-4) * (Tdsv[i] - 273.15)
    cpc_dsv = 1003.2 + 2.09 * (Tdsv[i] - 273.15)
    kc_dsv = 0.08 - (1e-4) * (Tdsv[i] - 273.15)

    # update wood and char density
    pw_dsv[i] = B_dsv[i]*rhow
    pc_dsv[i] = (C1_dsv[i]+C2_dsv[i])*rhow

    # update mass fraction vector
    Yw_dsv = pw_dsv[i] / (pw_dsv[i] + pc_dsv[i])
    cpbar_dsv = Yw_dsv*cpw_dsv + (1-Yw_dsv)*cpc_dsv
    kbar_dsv = Yw_dsv*kw_dsv + (1-Yw_dsv)*kc_dsv
    pbar_dsv = pw_dsv[i] + pc_dsv[i]
    Ys_dsv[i] = np.mean(B_dsv[i] + C1_dsv[i] + C2_dsv[i])

# Experimental Data from Sadhukhan 2009
# -----------------------------------------------------------------------------

t1, Tsph = np.loadtxt('sadhukhan2009/Fig2_Tcylinder.csv', delimiter=',', unpack=True)
t2, Msph = np.loadtxt('sadhukhan2009/Fig2_Mcylinder.csv', delimiter=',', unpack=True)

# Plot
#------------------------------------------------------------------------------

py.ion()
py.close('all')

def despine():
    ax = py.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    py.tick_params(axis='both', bottom='off', top='off', left='off', right='off')

py.figure(1)
py.plot(t, TdsvH[:, 0], 'b-', lw=2, label='$\Delta$H = -240 kJ/kg')
py.plot(t, Tdsv[:, 0], 'r-', lw=2, label='$\Delta$H = 0')
py.plot(t1, Tsph+273, 'go', mec='g', mew=2, label='experiment')
py.axhline(Tinf, c='k', ls='--')
py.ylim(ymin=Ti-20)
py.xlabel('Time (s)')
py.ylabel('Center Temperature (K)')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()

py.figure(2)
py.plot(t, Ys_dsvH, 'b-', lw=2, label='$\Delta$H = -240 kJ/kg')
py.plot(t, Ys_dsv, 'r-', lw=2, label='$\Delta$H = 0')
py.plot(t2, Msph, 'go', mec='g', mew=2, label='experiment')
py.ylim([0, 1.1])
py.xlabel('Time (s)')
py.ylabel('Residual Weight Fraction (-)')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()

# plot figures in black and white

py.figure(3)
py.plot(t, TdsvH[:, 0], c='k', ls='-', lw=2, label='$\Delta$H = -240 kJ/kg')
py.plot(t, Tdsv[:, 0], c='k', ls='--', lw=2, label='$\Delta$H = 0')
py.plot(t1, Tsph+273, c='k', marker='o', ls='', label='experiment')
py.axhline(Tinf, c='k', ls='--')
py.ylim(ymin=Ti-20)
py.xlabel('Time (s)')
py.ylabel('Center Temperature (K)')
py.legend(loc='best', numpoints=1, frameon=False)
despine()

py.figure(4)
py.plot(t, Ys_dsvH, c='k', ls='-', lw=2, label='$\Delta$H = -240 kJ/kg')
py.plot(t, Ys_dsv, c='k', ls='--', lw=2, label='$\Delta$H = 0')
py.plot(t2, Msph, c='k', marker='o', ls='', label='experiment')
py.ylim([0, 1.1])
py.xlabel('Time (s)')
py.ylabel('Residual Weight Fraction (-)')
py.legend(loc='best', numpoints=1, frameon=False)
despine()


