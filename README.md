# Code for Low-Order Particle Modeling Paper

Python code and COMSOL data for the following paper:

Gavin Wiggins, Peter Ciesielski, and Stuart Daw. "Low-Order Modeling of
Internal Heat Transfer in Biomass Particle Pyrolysis." Energy & Fuels, 2016,
30(6), pp. 4960-4969.  
[available online](http://pubs.acs.org/doi/abs/10.1021/acs.energyfuels.6b00554)

Requirements: Python 3, Matplotlib, NumPy, SciPy

## Data

Data from the 3-D Comsol simulations is available in the **comsol** folder.
Data from the Sadhukhan 2009 paper is available in the **sadhukhan2009**
folder.

## Functions

Module files prepended with **func** contain functions for modeling heat
conduction within a woody biomass particle. Details about each function are
available in the comment within each module file.

**funcHeatCond.py**  
Functions for 1-D transient heat conduction within a solid sphere, cylinder, or
slab shape. Each function returns an array of temperatures from the center to
surface of the particle at each time step. Properties such as thermal
conductivity and heat capacity can be constant or vary with temperature and
moisture content. Assumes convection at surface, symmetry at center, no
radiation, and constant particle size.

**funcKinetics.py**  
Kinetic reactions for gas, tar, and char yields from biomass pyrolysis.
Parameters for pre-factors and activation energies from Sadhukhan 2009 paper.

**funcOther.py**  
Various functions used to model 1-D biomass particle pyrolysis. Calculate the
shell volumes that comprise a solid sphere. Calculate the volume average
temperature of the entire sphere. Calculate the Sauter diameter of a shape.
Calculate the dimensionless Biot and pyrolysis numbers.

**funcRoots.py , funcTheta.py , funcZeta.py**  
Functions for solving the 1-D analytical solution of the transient heat
conduction equation.

## Models

Various models were created to investigate internal heat transfer within a
solid woody biomass particle. See the comments in each file for detailed
documentation, alternatively an overview is provided below.

**oak-200-20000.py**  
Compare volume average temperature profiles from 1-D model and 3-D Comsol
simulation of white oak particles with Feret diameters DF = 200 um to 20 mm.
Surface area to volume diameter, Dsv, is used for the 1-D model.

**oak-200.py , oak-20000.py**  
Compare temperature profiles of 1-D and 3-D models for DF = 200 um and 20 mm
dry white oak particle. Heat capacity as function of temperature and constant
thermal conductivity. Different equivalent spherical diameters and
characteristic lengths implemented with 1-D model.

**oak-bipy.py**  
Compare Biot and pyrolysis numbers for dry white oak particles. Thermal
properties evaluated at 773 K and kinetic rate constant from Sadhukhan 2009
paper.

**oak-diff.py**
Calculate and compare difference between 3-D surface and center temperature
profiles.

**pine-200.py , pine-20000.py**  
Compare temperature profiles of 1-D and 3-D models for DF = 200 um and 20 mm
dry loblolly pine particle. Heat capacity as function of temperature and
constant thermal conductivity. Different equivalent spherical diameters and
characteristic lengths implemented with 1-D model.

**pine-5400.py**  
Surface, volume, and center temperatures from 1-D and 3-D models for DF = 5.4
mm dry loblolly pine particle. Heat capacity as function of temperature and
constant thermal conductivity. Dsv used as equivalent spherical diameter for
1-D model.

**sadhukhan2009.py**  
Compare 1-D transient heat conduction model to Sadhukhan 2009 Figure 2
cylinder.

**sphere-ana-num.py**  
Compare 1-D analytical sphere solution to 1-D numerical and 3-D Comsol
solutions for transient heat conduction in solid sphere with constant k and Cp.

**sphere-cube.py**  
Compare temperature profiles of 3-D cube and 3-D sphere in Comsol to 1-D cube
and 1-D sphere model. Compare 3-D cube and 3-D sphere Comsol temperature
profiles. Sphere and cube are volume equivalent where sphere diameter is 1 mm
and cube side is 0.806 mm.

**sphere-cyl-slab.py**  
Compare temperature profiles of 1-D solid sphere, cylinder, and cube shapes
that are volume equivalent. Note that due to surface area, sphere heats the
slowest compared to the cylinder and cube shapes.

## License

Files in this repository are available under the MIT license. See the LICENSE
file for more info.

