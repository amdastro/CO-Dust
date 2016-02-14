import physical as ph
import numpy as np

'''----------RUNTIME PARAMETERS--------'''



# Set the default timestep, and the total simulation time
#dt_init = 1e3
#tmin = 1e4
#tmax = 1e8
tmin = 3.*ph.daytosec
tmax = 50.*ph.daytosec
dt_init = 1e1
#integrate after this time (later set it to integrate at given some physical conditions)
t_integrate = tmin

# Limit the percent of CO that can form in each timestep:
percentCO = 1e-3
# The amount of CO above which we apply the percent criteria
CO_min = 1e-2

# The total mass fraction adds to 1
# X_i = A_i * Y_i
# Evolve things in terms of number fraction (Y is number of molecules)
# but view them in mass fraction
X_C_tot = 0.05
X_O_tot = 0.1
Y_C_tot = X_C_tot / ph.A_C
Y_O_tot = X_O_tot / ph.A_O

directory = 'dt%.0e_min%.0e_per%.0e_XC%.2f_XO%.2f_PR'%(dt_init,percentCO,CO_min,X_C_tot,X_O_tot)


# The total mass of the cold shell - given by the gamma ray luminosity:
# (sets the amount of nonthermal electrons)
M_cs = 1e-4*ph.Msun 

''' The initial thermodynamics of the cold shell '''
# The ejecta velocity: 
#v_ej = 1e8 
# The initial temperature, assuming it has already cooled efficiently
#T_cs_init = 1e4
# The initial radius of shell
#R_cs_init = 5e13
# Initial thickness:
#delta_init = R_cs_init/4000.
# consequent initial density:
#n_init = M_cs/(4.*np.pi * ph.mp * R_cs_init**2 * delta_init)
# Adiabatic index
#gamma = 1.3


''' The initial thermodynamics of Pontefract & Rawlings'''
# The initial temperature
T_cs_init = 6000.
# initial density:
n_init = 1.e15
# photospheric radius
R_cs_init = 1.44e12
# ejecta velocity
v_ej = 3.21e12/ph.daytosec


'''To read trajectory from textfile
Read in now and interpolate from arrays 
in trajectory.py'''
#t_array, n_array, T_array = np.genfromtxt("traj_shocklum.dat",unpack=True)


# Gamma rays:
# The observed gamma ray luminosity:
L_gamma = 3e35 
epsilon = 1.
# The energy per particle produced
W_d = 2.4e-10


# floor value (so things don't go to zero)
floor = 1e-100