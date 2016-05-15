import physical as ph
import numpy as np

'''----------RUNTIME PARAMETERS--------'''



# Set the default timestep
dt_init = 1e3

# Limit the percent change of Y in each timestep:
percent = 1e-2
# The minimum Y above which we apply the percent criteria
Y_min = 1e-5

# The total mass fraction adds to 1
# X_i = A_i * Y_i
# Evolve things in terms of number fraction (Y is number of molecules)
#X_C_tot  = 0.05
#X_O_tot  = 0.1
#X_Si_tot = 6e-4
#X_Mg_tot = 5e-4
#X_other  = 1. - X_C_tot - X_O_tot - X_Si_tot - X_Mg_tot
'''
X_H_tot  = 0.75
#X_He_tot = 0.394
X_C_tot  = 0.05
#X_N_tot  = 0.125
X_O_tot  = 0.1
#X_Ne_tot = 2.3e-3
X_Mg_tot = 0.005
X_Si_tot = 0.005
'''

# Morisset & Pequignot 1996 abundances
X_H_tot  = 0.372
#X_He_tot = 0.394
X_C_tot  = 8.0e-3
#X_N_tot  = 0.125
X_O_tot  = 9.5e-2 
#X_Ne_tot = 2.3e-3
X_Mg_tot = 6.7e-4
X_Si_tot = 7.7e-4

#X_S_tot  = 5.4e-4
#X_Cl_tot = 2.4e-5
#X_Ar_tot = 1.6e-4
#X_Ca_tot = 1.9e-4
#X_Fe_tot = 1.6e-3

# This adds to 1.00023 ...

Y_H_tot  =  X_H_tot  / ph.A_H
#Y_He_tot =  X_He_tot / ph.A_He
Y_C_tot  =  X_C_tot  / ph.A_C
#Y_N_tot  =  X_N_tot  / ph.A_N
Y_O_tot  =  X_O_tot  / ph.A_O
#Y_Ne_tot =  X_Ne_tot / ph.A_Ne
Y_Mg_tot =  X_Mg_tot / ph.A_Mg
Y_Si_tot =  X_Si_tot / ph.A_Si
#Y_S_tot  =  X_S_tot  / ph.A_S
#Y_Cl_tot =  X_Cl_tot / ph.A_Cl
#Y_Ar_tot =  X_Ar_tot / ph.A_Ar
#Y_Ca_tot =  X_Ca_tot / ph.A_Ca
#Y_Fe_tot =  X_Fe_tot / ph.A_Fe

# Mass fractions of species not included in chemical reactions
#X_other = X_H_tot + X_He_tot + X_N_tot + X_S_tot + X_Cl_tot + X_Ar_tot + X_Ca_tot + X_Fe_tot
#Y_other = Y_H_tot + Y_He_tot + Y_N_tot + Y_S_tot + Y_Cl_tot + Y_Ar_tot + Y_Ca_tot + Y_Fe_tot
X_other = 1. - X_C_tot - X_O_tot - X_Mg_tot - X_Si_tot - X_H_tot
Y_other = X_other # NOT REALLY unless it's all Hydrogen! which it's totally not. 

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
in integrate.py'''
t_array, n_array, T_array, L_array = np.genfromtxt("traj_shocklum.dat",unpack=True)
n_init = n_array[0]
T_cs_init = T_array[0]
tmin = t_array[0]
tmax = t_array[-1]


# when to begin integrating instead of assuming chemical equilibrium
t_integrate = tmin

# Gamma rays:
# The observed gamma ray luminosity - this now declines in time, read from table
L_shock_init = L_array[0]
epsilon = 1e-3
# The energy per particle produced
# Currently both are 125 eV
W_d_CO  = 2.4e-10
W_d_SiO = 2.4e-10

directory = 'MP/dt1e3_per1e-2/shocklum_e%.0e'%(epsilon)

# length of dust arrays
arraylen = 200000
# floor value (so things don't go to zero)
floor = 1e-200