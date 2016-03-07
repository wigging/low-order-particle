"""
Compare Biot and pyrolysis numbers for dry white oak particles. Thermal properties
evaluated at 773 K and kinetic rate constant from Sadhukhan 2009 paper.
"""
import numpy as np
import matplotlib.pyplot as py
from funcHeatCond import hc2
from funcOther import vol, Tvol, biot, py1, py2

# Parameters
# -----------------------------------------------------------------------------

Gb = 0.72       # basic specific gravity, Wood Handbook Table 4-7, (-)
rho = Gb*1000   # density, kg/m^3
k = 0.16        # thermal conductivity, W/mK
cp = 3092       # heat capacity, J/kgK
h = 350         # heat transfer coefficient, W/m^2*K
Ti = 293        # initial particle temp, K
Tinf = 773      # ambient temp, K

# Kinetic Parameters from Sadhukhan 2009
# -----------------------------------------------------------------------------

# A as pre-factor (1/s) and E as activation energy (kJ/mol) from Sadhukhan 2009
A1 = 168.4;     E1 = 51.965;    # biomass -> (vol+gas)
A2 = 13.2;      E2 = 45.960;    # biomass -> char

# evaluate reaction rate constant for wood conversion
R = 0.008314                        # universal gas constant, kJ/mol*K
K1 = A1 * np.exp(-E1 / (R*Tinf))    # biomass -> (vol+gas)
K2 = A2 * np.exp(-E2 / (R*Tinf))    # biomass -> char
kr = K1 + K2                        # rate constant for wood conversion, 1/s

# Particle Size of DF = 200 um
# -----------------------------------------------------------------------------

# particle size data from Comsol for DF = 200 um
As200 = 5.355e-8   # surface area of Comsol particle, m^2
v200 = 8.895e-13   # volume of Comsol particle, m^3
H200 = 1.8914e-4   # height of Comsol particle, m
W200 = 6.9766e-5   # width of Comsol particle, m
L200 = 8.8923e-5   # length of Comsol particle, m

# equivalent spherical diameters and characteristic lengths for DF = 200um
ds200 = (As200/np.pi)**(1/2)    # surface area equivalent sphere diameter, m
dv200 = (6/np.pi*v200)**(1/3)   # volume equivalent sphere diameter, m
dsv200 = (dv200**3)/(ds200**2)  # surface volume equivalent sphere diameter (Sauter), m
dc200 = v200/As200              # characteristic length, m
dh200 = H200                    # height as diameter, m
dw200 = W200                    # width as diameter, m
dl200 = L200                    # length as diameter, m

# Biot numbers for DF = 200 um, Bi (-)
Bi_ds200 = biot(h, ds200/2, k)
Bi_dv200 = biot(h, dv200/2, k)
Bi_dsv200 = biot(h, dsv200/2, k)
Bi_dc200 = biot(h, dc200/2, k)
Bi_dh200 = biot(h, dh200/2, k)
Bi_dw200 = biot(h, dw200/2, k)
Bi_dl200 = biot(h, dl200/2, k)

# Pyrolysis numbers as Py I and Py II for DF = 200 um, Py (-)
Py1_ds200 = py1(k, kr, rho, cp, ds200/2)
Py2_ds200 = py2(h, kr, rho, cp, ds200/2)

Py1_dv200 = py1(k, kr, rho, cp, dv200/2)
Py2_dv200 = py2(h, kr, rho, cp, dv200/2)

Py1_dsv200 = py1(k, kr, rho, cp, dsv200/2)
Py2_dsv200 = py2(h, kr, rho, cp, dsv200/2)

Py1_dc200 = py1(k, kr, rho, cp, dc200/2)
Py2_dc200 = py2(h, kr, rho, cp, dc200/2)

Py1_dh200 = py1(k, kr, rho, cp, dh200/2)
Py2_dh200 = py2(h, kr, rho, cp, dh200/2)

Py1_dw200 = py1(k, kr, rho, cp, dw200/2)
Py2_dw200 = py2(h, kr, rho, cp, dw200/2)

Py1_dl200 = py1(k, kr, rho, cp, dl200/2)
Py2_dl200 = py2(h, kr, rho, cp, dl200/2)

# Particle Size of DF = 400 um
# -----------------------------------------------------------------------------

# particle size data from Comsol for DF = 400 um
As400 = 1.879e-7   # surface area of Comsol particle, m^2
v400 = 5.553e-12   # volume of Comsol particle, m^3

# equivalent spherical diameters for DF = 400 um
ds400 = (As400/np.pi)**(1/2)    # surface area equivalent sphere diameter, m
dv400 = (6/np.pi*v400)**(1/3)   # volume equivalent sphere diameter, m
dsv400 = (dv400**3)/(ds400**2)  # surface volume equivalent sphere diameter (Sauter), m

# Biot and Pyrolysis Numbers for DF = 400 um
Bi_dsv400 = biot(h, dsv400/2, k)
Py1_dsv400 = py1(k, kr, rho, cp, dsv400/2)
Py2_dsv400 = py2(h, kr, rho, cp, dsv400/2)

# Particle Size of DF = 700 um
# -----------------------------------------------------------------------------

# particle size data from Comsol for DF = 700 um
As700 = 4.836e-7   # surface area of Comsol particle, m^2
v700 = 2.11e-11    # volume of Comsol particle, m^3

# equivalent spherical diameters for DF = 700 um
ds700 = (As700/np.pi)**(1/2)    # surface area equivalent sphere diameter, m
dv700 = (6/np.pi*v700)**(1/3)   # volume equivalent sphere diameter, m
dsv700 = (dv700**3)/(ds700**2)  # surface volume equivalent sphere diameter (Sauter), m

# Biot and Pyrolysis Numbers for DF = 700 um
Bi_dsv700 = biot(h, dsv700/2, k)
Py1_dsv700 = py1(k, kr, rho, cp, dsv700/2)
Py2_dsv700 = py2(h, kr, rho, cp, dsv700/2)

# Particle Size of DF = 1400 um
# -----------------------------------------------------------------------------

# particle size data from Comsol for DF = 1400 um
As1400 = 1.394e-6   # surface area of Comsol particle, m^2
v1400 = 8.442e-11   # volume of Comsol particle, m^3

# equivalent spherical diameters for DF = 1400 um
ds1400 = (As1400/np.pi)**(1/2)    # surface area equivalent sphere diameter, m
dv1400 = (6/np.pi*v1400)**(1/3)   # volume equivalent sphere diameter, m
dsv1400 = (dv1400**3)/(ds1400**2)  # surface volume equivalent sphere diameter (Sauter), m

# Biot and Pyrolysis Numbers for DF = 1400 um
Bi_dsv1400 = biot(h, dsv1400/2, k)
Py1_dsv1400 = py1(k, kr, rho, cp, dsv1400/2)
Py2_dsv1400 = py2(h, kr, rho, cp, dsv1400/2)

# Particle Size of DF = 2800 um
# -----------------------------------------------------------------------------

# particle size data from Comsol for DF = 2800 um
As2800 = 4.614e-6   # surface area of Comsol particle, m^2
v2800 = 4.011e-10   # volume of Comsol particle, m^3

# equivalent spherical diameters for DF = 2800 um
ds2800 = (As2800/np.pi)**(1/2)    # surface area equivalent sphere diameter, m
dv2800 = (6/np.pi*v2800)**(1/3)   # volume equivalent sphere diameter, m
dsv2800 = (dv2800**3)/(ds2800**2)  # surface volume equivalent sphere diameter (Sauter), m

# Biot and Pyrolysis Numbers for DF = 2800 um
Bi_dsv2800 = biot(h, dsv2800/2, k)
Py1_dsv2800 = py1(k, kr, rho, cp, dsv2800/2)
Py2_dsv2800 = py2(h, kr, rho, cp, dsv2800/2)

# Particle Size of DF = 5400 um
# -----------------------------------------------------------------------------

# particle size data from Comsol for DF = 5400 um
As5400 = 1.716e-5   # surface area of Comsol particle, m^2
v5400 = 2.877e-9    # volume of Comsol particle, m^3

# equivalent spherical diameters for DF = 5400 um
ds5400 = (As5400/np.pi)**(1/2)    # surface area equivalent sphere diameter, m
dv5400 = (6/np.pi*v5400)**(1/3)   # volume equivalent sphere diameter, m
dsv5400 = (dv5400**3)/(ds5400**2)  # surface volume equivalent sphere diameter (Sauter), m

# Biot and Pyrolysis Numbers for DF = 5400 um
Bi_dsv5400 = biot(h, dsv5400/2, k)
Py1_dsv5400 = py1(k, kr, rho, cp, dsv5400/2)
Py2_dsv5400 = py2(h, kr, rho, cp, dsv5400/2)

# Particle Size of DF = 10000 um
# -----------------------------------------------------------------------------

# particle size data from Comsol for DF = 10000 um
As10000 = 5.885e-5   # surface area of Comsol particle, m^2
v10000 = 1.827e-8    # volume of Comsol particle, m^3

# equivalent spherical diameters for DF = 10000 um
ds10000 = (As10000/np.pi)**(1/2)    # surface area equivalent sphere diameter, m
dv10000 = (6/np.pi*v10000)**(1/3)   # volume equivalent sphere diameter, m
dsv10000 = (dv10000**3)/(ds10000**2)  # surface volume equivalent sphere diameter (Sauter), m

# Biot and Pyrolysis Numbers for DF = 10000 um
Bi_dsv10000 = biot(h, dsv10000/2, k)
Py1_dsv10000 = py1(k, kr, rho, cp, dsv10000/2)
Py2_dsv10000 = py2(h, kr, rho, cp, dsv10000/2)

# Particle Size of DF = 20000 um
# -----------------------------------------------------------------------------

# particle size data from Comsol for DF = 20000 um
As20000 = 2.354e-4   # surface area of Comsol particle, m^2
v20000 = 1.462e-7    # volume of Comsol particle, m^3
H20000 = 0.02027     # height of Comsol particle, m
W20000 = 0.001877    # width of Comsol particle, m
L20000 = 0.005069    # length of Comsol particle, m

# equivalent spherical diameters and characteristic lengths for DF = 200um
ds20000 = (As20000/np.pi)**(1/2)        # surface area equivalent sphere diameter, m
dv20000 = (6/np.pi*v20000)**(1/3)       # volume equivalent sphere diameter, m
dsv20000 = (dv20000**3)/(ds20000**2)    # surface volume equivalent sphere diameter (Sauter), m
dc20000 = v20000/As20000                # characteristic length, m
dh20000 = H20000                        # height as diameter, m
dw20000 = W20000                        # width as diameter, m
dl20000 = L20000                        # length as diameter, m

# Biot numbers for DF = 20 mm, Bi (-)
Bi_ds20000 = biot(h, ds20000/2, k)
Bi_dv20000 = biot(h, dv20000/2, k)
Bi_dsv20000 = biot(h, dsv20000/2, k)
Bi_dc20000 = biot(h, dc20000/2, k)
Bi_dh20000 = biot(h, dh20000/2, k)
Bi_dw20000 = biot(h, dw20000/2, k)
Bi_dl20000 = biot(h, dl20000/2, k)

# Pyrolysis numbers as Py I and Py II for DF = 20 mm, Py (-)
Py1_ds20000 = py1(k, kr, rho, cp, ds20000/2)
Py2_ds20000 = py2(h, kr, rho, cp, ds20000/2)

Py1_dv20000 = py1(k, kr, rho, cp, dv20000/2)
Py2_dv20000 = py2(h, kr, rho, cp, dv20000/2)

Py1_dsv20000 = py1(k, kr, rho, cp, dsv20000/2)
Py2_dsv20000 = py2(h, kr, rho, cp, dsv20000/2)

Py1_dc20000 = py1(k, kr, rho, cp, dc20000/2)
Py2_dc20000 = py2(h, kr, rho, cp, dc20000/2)

Py1_dh20000 = py1(k, kr, rho, cp, dh20000/2)
Py2_dh20000 = py2(h, kr, rho, cp, dh20000/2)

Py1_dw20000 = py1(k, kr, rho, cp, dw20000/2)
Py2_dw20000 = py2(h, kr, rho, cp, dw20000/2)

Py1_dl20000 = py1(k, kr, rho, cp, dl20000/2)
Py2_dl20000 = py2(h, kr, rho, cp, dl20000/2)

# Print Biot Numbers
# -----------------------------------------------------------------------------

# Bi < 1 with Py II, Biot number less than 1 compare with Py II
# Bi > 1 with Py I, Biot numbers less greater than 1 compare with Py I

# Biot numbers for DF = 200 um, Bi (-)
print('DF = 200 um')
print('Bi ds200 =', Bi_ds200)
print('Bi dv200 =', Bi_dv200)
print('Bi dsv200 =', Bi_dsv200)
print('Bi dc200 =', Bi_dc200)
print('Bi dh200 =', Bi_dh200)
print('Bi dw200 =', Bi_dw200)
print('Bi dl200 =', Bi_dl200, '\n')

# Biot numbers for Dsv = 200 um to 20 mm
print('Dsv = 200 um to 20 mm')
print('Bi dsv200 =', Bi_dsv200)
print('Bi dsv400 =', Bi_dsv400)
print('Bi dsv700 =', Bi_dsv700)
print('Bi dsv1400 =', Bi_dsv2800)
print('Bi dsv2800 =', Bi_dsv2800)
print('Bi dsv5400 =', Bi_dsv5400)
print('Bi dsv10000 =', Bi_dsv10000)
print('Bi dsv20000 =', Bi_dsv20000, '\n')

# Biot numbers for DF = 20 mm, Bi (-)
print('DF = 20 mm')
print('Bi ds20000 =', Bi_ds20000)
print('Bi dv20000 =', Bi_dv20000)
print('Bi dsv20000 =', Bi_dsv20000)
print('Bi dc20000 =', Bi_dc20000)
print('Bi dh20000 =', Bi_dh20000)
print('Bi dw20000 =', Bi_dw20000)
print('Bi dl20000 =', Bi_dl20000)

# Plot Results
# -----------------------------------------------------------------------------

py.ion()
py.close('all')

def despine():
    ax = py.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    py.tick_params(axis='both', bottom='on', top='off', left='on', right='off')
    py.minorticks_off()

def addtext():
    ax = py.gca()
    ax.text(0.2, 0.93, 'Py II\nkinetics, isothermal', ha='center', transform=ax.transAxes)
    ax.text(0.8, 0.93, 'Py I\nkinetics', ha='center', transform=ax.transAxes)
    ax.text(0.2, 0.03, 'Py II\nconvection', ha='center', transform=ax.transAxes)
    ax.text(0.8, 0.03, 'Py I\nconduction', ha='center', transform=ax.transAxes)

py.figure(1, figsize=(20, 6))

py.subplot(1, 2, 1)
py.subplots_adjust(wspace=0.1)
py.plot(Bi_ds20000, Py1_ds20000, 'o', ms=10, mec='none', label='$\mathregular{D_S}$')
py.plot(Bi_dv20000, Py1_dv20000, 'o', ms=10, mec='none', label='$\mathregular{D_V}$')
py.plot(Bi_dsv20000, Py1_dsv20000, 'o', ms=10, mec='none', label='$\mathregular{D_{SV}}$')
py.plot(Bi_dc20000, Py2_dc20000, 'o', ms=10, mec='none', label='$\mathregular{D_{CH}}$')
py.plot(Bi_dh20000, Py1_dh20000, 'o', ms=10, mec='none', label='$\mathregular{D_H}$')
py.plot(Bi_dw20000, Py1_dw20000, 'o', ms=10, mec='none', label='$\mathregular{D_W}$')
py.plot(Bi_dl20000, Py1_dl20000, 'o', ms=10, mec='none', label='$\mathregular{D_L}$')
py.plot(Bi_ds200, Py2_ds200, '^', ms=10, mec='none')
py.plot(Bi_dv200, Py2_dv200, '^', ms=10, mec='none')
py.plot(Bi_dsv200, Py2_dsv200, '^', ms=10, mec='none')
py.plot(Bi_dc200, Py2_dc200, '^', ms=10, mec='none')
py.plot(Bi_dh200, Py2_dh200, '^', ms=10, mec='none')
py.plot(Bi_dw200, Py2_dw200, '^', ms=10, mec='none')
py.plot(Bi_dl200, Py2_dl200, '^', ms=10, mec='none')
py.plot([], [], '^', ms=10, mew=1, mfc='none', label='0.2 mm')
py.plot([], [], 'o', ms=10, mew=1, mfc='none', label='20 mm')
py.xlabel('Biot Number, Bi (-)')
py.ylabel('Pyrolysis Number, Py (-)')
py.title('DF =  200 um (triangle) and 20 mm (circle)')
py.axvline(1, c='k', ls='-.')
py.axvspan(10**-1, 10**1, color='0.9')
py.axhline(1, c='k', ls='-.')
py.axhspan(10**-1, 10**1, color='0.9')
py.xlim(10**-4, 10**4)
py.ylim(10**-4, 10**4)
py.xscale('log')
py.yscale('log')
py.legend(loc='center right', numpoints=1, frameon=False)
despine()
addtext()

py.subplot(1, 2, 2)
py.plot(Bi_dsv200, Py2_dsv200, 'o', ms=10, mec='none', label='0.2 mm')
py.plot(Bi_dsv400, Py2_dsv400, 'o', ms=10, mec='none', label='0.4 mm')
py.plot(Bi_dsv700, Py2_dsv700, 'o', ms=10, mec='none', label='0.7 mm')
py.plot(Bi_dsv1400, Py2_dsv1400, 'o', ms=10, mec='none', label='1.4 mm')
py.plot(Bi_dsv2800, Py2_dsv2800, 'o', ms=10, mec='none', label='2.8 mm')
py.plot(Bi_dsv5400, Py1_dsv5400, 'o', ms=10, mec='none', label='5.4 mm')
py.plot(Bi_dsv10000, Py1_dsv10000, 'o', ms=10, mec='none', label='10 mm')
py.plot(Bi_dsv20000, Py1_dsv20000, 'o', ms=10, mec='none', label='20 mm')
py.xlabel('Biot Number, Bi (-)')
py.title('Dsv =  200 um - 20 mm')
py.axvline(1, c='k', ls='-.')
py.axvspan(10**-1, 10**1, color='0.9')
py.axhline(1, c='k', ls='-.')
py.axhspan(10**-1, 10**1, color='0.9')
py.xlim(10**-4, 10**4)
py.ylim(10**-4, 10**4)
py.xscale('log')
py.yscale('log')
py.legend(loc='center right', numpoints=1, frameon=False, borderpad=0)
ax = py.gca()
ax.axes.get_yaxis().set_visible(False)
despine()
addtext()
