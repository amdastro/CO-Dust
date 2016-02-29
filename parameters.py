import physical as ph
import numpy as np

'''----------RUNTIME PARAMETERS--------'''



# Set the default timestep
#dt_init = 1e3
#tmin = 1e4
#tmax = 1e8
dt_init = 1e3

# Limit the percent change of Y in each timestep:
percent = 1e-3
# The amount of CO above which we apply the percent criteria
Y_min = 1e-2

# The total mass fraction adds to 1
# X_i = A_i * Y_i
# Evolve things in terms of number fraction (Y is number of molecules)
# but view them in mass fraction
X_C_tot = 0.05
X_O_tot = 0.1
Y_C_tot = X_C_tot / ph.A_C
Y_O_tot = X_O_tot / ph.A_O

# The total mass of the cold shell - given by the gamma ray luminosity:
# (sets the amount of nonthermal electrons)
M_cs = 1e-4*ph.Msun 

''' The initial thermodynamics of Pontefract & Rawlings'''
# The initial temperature
#T_cs_init = 6000.
# initial density:
#n_init = 1.e15
# photospheric radius
#R_cs_init = 1.44e12
# ejecta velocity
#v_ej = 3.21e12/ph.daytosec
#tmin = 3.*ph.daytosec
#tmax = 50.*ph.daytosec

'''To read trajectory from textfile
Read in now and interpolate from arrays 
in trajectory.py'''
t_array, n_array, T_array, L_array = np.genfromtxt("traj_shocklum2.dat",unpack=True)
n_init = n_array[0]
T_cs_init = T_array[0]
tmin = t_array[0]
tmax = t_array[-1]


# when to begin integrating instead of assuming chemical equilibrium
t_integrate = tmin

# Gamma rays:
# The observed gamma ray luminosity - this now declines in time, read from table
L_shock_init = L_array[0]
epsilon = 0.001
# The energy per particle produced
W_d = 2.4e-10

directory = 'shocklum_dt%.0e_XC%.2f_XO%.2f_e%.0e'%(dt_init,X_C_tot,X_O_tot,epsilon)

# length of dust arrays
arraylen = 200000
# floor value (so things don't go to zero)
floor = 1e-200