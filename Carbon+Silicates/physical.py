import numpy as np


mp    = 1.67e-24
kB    = 1.38e-16
h     = 6.626e-27
Msun  = 1.9e33

muCO  = 1.15e-23 # reduced mass of C and O
muSiO = 1.70e-23
B_CO  = 1.78e-11
B_SiO = 1.33e-11

# carbon dust:
mC = 12*1.67e-24 	# mass of a carbon atom
cshape = (36*np.pi)**(1./3.) # Shape factor for sphere. This is the ratio of the Survace to the Volume^(2/3). 
rhoC = 2.23 	# Density of solid carbon
sigmaC = 1500. # Surface energy of solid carbon
vC = mC/rhoC  # Volume occupied by a carbon atom in the bulk of the solid phase

# silicate dust:
'''UPDATE THESE!!'''
mMg2SiO4 = 12*1.67e-24 	# mass of Mg2SiO4 molecule
rhoMg2SiO4 = 2.23 	# Density of solid carbon
sigmaMg2SiO4 = 1500. # Surface energy of solid carbon
vMg2SiO4 = mMg2SiO4/rhoMg2SiO4  # Volume occupied by a carbon atom in the bulk of the solid phase

# atomic numbers for mass fraction: 
A_H  = 1.01
A_He = 4.0
A_C  = 12.011
A_N  = 14.01
A_O  = 15.999
A_Ne = 20.18
A_Mg = 24.205
A_Si = 28.086
A_S  = 32.07
A_Cl = 35.45
A_Ar = 39.95
A_Ca = 40.08
A_Fe = 55.85
A_CO = A_C + A_O 
A_SiO = A_Si + A_O
A_Mg2SiO4 = 2.*A_Mg + A_SiO + 3.*A_O

# conversion parameters
daytosec = 86400.

# parameters for Gibbs free energy of Mg2SiO4
# From Nozawa et al.
GibbsA = 18.62 * 1e4 #K
GibbsB = 52.4336