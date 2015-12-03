import physical as ph
import numpy as np

'''----------RUNTIME PARAMETERS--------'''


# Set the default timestep, and the total simulation time
dt_init = 1e3
tmin = 1e4
tmax = 1e6
#integrate after this time (later set it to integrate at given some physical conditions)
t_integrate = 1.e5

# Limit the percent of CO that can form in each timestep:
percentCO = 1e-3
# The amount of CO above which we apply the percent criteria
CO_amount = 1e-3

directory = 'tint_%.0e_dt_%.0e_percent%.0e_amt%.0e'%(t_integrate,dt_init,percentCO,CO_amount)


''' The initial thermodynamics of the cold shell '''
# The total mass of the cold shell - given by the gamma ray luminosity:
M_cs = 2e-4*ph.Msun 
# The ejecta velocity: 
v_ej = 1e8 
# The initial temperature, assuming it has already cooled efficiently
T_cs_init = 1e4
# The initial radius of shell
R_cs_init = 1e14
# Initial thickness:
delta_init = R_cs_init/4000.
# consequent initial density:
n_init = M_cs/(4.*np.pi * ph.mp * R_cs_init**2 * delta_init)

# Adiabatic index
gamma = 1.3


# The total number fraction of C and O atoms:
Y_C_tot = 0.4
Y_O_tot = 0.6


# Gamma rays:
# The observed gamma ray luminosity:
L_gamma = 3e35 
epsilon = 1.
# The energy per particle produced
W_d = 2.4e-10