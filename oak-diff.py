"""
Calculate and compare difference between 3-D surface and center temperature profiles.
"""

import numpy as np
import matplotlib.pyplot as py

# Data from Comsol 3-D Oak Particle Simulation
# -----------------------------------------------------------------------------

f200 = 'comsol/200tempsOak.txt'
t200, _, _, Tc200, _, _, Tsa200 = np.loadtxt(f200, skiprows=5, unpack=True)
diff200 = abs(Tsa200 - Tc200)
max200 = np.max(diff200)

f400 = 'comsol/400tempsOak.txt'
t400, _, _, Tc400, _, _, Tsa400 = np.loadtxt(f400, skiprows=5, unpack=True)
diff400 = abs(Tsa400 - Tc400)
max400 = np.max(diff400)

f700 = 'comsol/700tempsOak.txt'
t700, _, _, Tc700, _, _, Tsa700 = np.loadtxt(f700, skiprows=5, unpack=True)
diff700 = abs(Tsa700 - Tc700)
max700 = np.max(diff700)

f1400 = 'comsol/1400tempsOak.txt'
t1400, _, _, Tc1400, _, _, Tsa1400 = np.loadtxt(f1400, skiprows=5, unpack=True)
diff1400 = abs(Tsa1400 - Tc1400)
max1400 = np.max(diff1400)

f2800 = 'comsol/2800tempsOak.txt'
t2800, _, _, Tc2800, _, _, Tsa2800 = np.loadtxt(f2800, skiprows=5, unpack=True)
diff2800 = abs(Tsa2800 - Tc2800)
max2800 = np.max(diff2800)

f5400 = 'comsol/5400tempsOak.txt'
t5400, _, _, Tc5400, _, _, Tsa5400 = np.loadtxt(f5400, skiprows=5, unpack=True)
diff5400 = abs(Tsa5400 - Tc5400)
max5400 = np.max(diff5400)

f10000 = 'comsol/10000tempsOak.txt'
t10000, _, _, Tc10000, _, _, Tsa10000 = np.loadtxt(f10000, skiprows=5, unpack=True)
diff10000 = abs(Tsa10000 - Tc10000)
max10000 = np.max(diff10000)

f20000 = 'comsol/20000tempsOak.txt'
t20000, _, _, Tc20000, _, _, Tsa20000 = np.loadtxt(f20000, skiprows=5, unpack=True)
diff20000 = abs(Tsa20000 - Tc20000)
max20000 = np.max(diff20000)

# bar chart
maxDiff = [max200, max400, max700, max1400, max2800, max5400, max10000, max20000]
xlabels = ('200um', '400um', '700um', '1.4mm', '2.8mm', '5.4mm', '10mm', '20mm')
xloc = range(len(xlabels))

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
py.plot(t20000, diff20000, lw=2, label='20 mm')
py.plot(t10000, diff10000, lw=2, label='10 mm')
py.plot(t5400, diff5400, lw=2, label='5.4 mm')
py.plot(t2800, diff2800, lw=2, label='2.8 mm')
py.plot(t1400, diff1400, lw=2, label='1.4 mm')
py.plot(t700, diff700, lw=2, label='700 um')
py.plot(t400, diff400, lw=2, label='400 um')
py.plot(t200, diff200, lw=2, label='200 um')
py.xlim(0, 40)
py.ylabel('Temperature Difference (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1, frameon=False)
py.grid()
despine()

py.figure(2)
py.bar(xloc, maxDiff, align='center')
py.xticks(xloc, xlabels)
py.ylabel('Max Temperature Difference (K)')
py.xlabel('Feret Diameter')
py.grid()
despine()
