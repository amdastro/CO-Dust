import numpy as np


mp = 1.67e-24
kB = 1.38e-16
h = 6.626e-27
mu = 1.e-23 # reduced mass of C and O
B_CO = 1.78e-11
Msun = 1.9e33

# dust:
mC = 12*1.67e-24 	# mass of a carbon atom
cshape = (36*np.pi)**(1./3.) # Shape factor for sphere. This is the ratio of the Survace to the Volume^(2/3). 
rhoC = 2.23 	# Density of solid carbon
sigma = 1500. # Surface energy of solid carbon
vC = mC/rhoC  # Volume occupied by a carbon atom in the bulk of the solid phase

# atomic numbers for mass fraction: 
A_C = 12.
A_O = 15.999
A_CO =A_C+A_O # just about?

# conversion parameters
daytosec = 86400.
