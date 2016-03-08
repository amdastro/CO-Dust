import numpy as np
import pylab as plt
import physical as ph
import reaction as re
import parameters as par
import abundances as ab
import trajectory as traj
from timestep import choosedt
from scipy.interpolate import interp1d
import os
import sys

print 'RUNNING ', par.directory

#----------- Initial fluid parameters ------------
T_cs = par.T_cs_init
n = par.n_init
# interpolated arrays
n_interp = interp1d(par.t_array,par.n_array,kind='linear')
T_interp = interp1d(par.t_array,par.T_array,kind='linear')
L_interp = interp1d(par.t_array,par.L_array,kind='linear')

#------------ Initial reaction rates ------------
# L_gamma is set by shock power L_shock which decreases in time
# this is self consistently calcualted with the trajectory
K_ra = re.formation(T_cs)
K_th = re.thermal(T_cs)
# K_rd is the sum of both destruction rates. 
# f_CO is the fraction of electrons that impact CO
#f_CO = Y_CO/(Y_CO+Y_C_free+Y_O_free+Y_D+par.Y_other)
# ? but haven't set other abundances yet?
f_CO=0.
K_nth = par.L_shock_init*par.epsilon*ph.mp*f_CO/par.M_cs/par.W_d
K_rd = K_th + K_nth

#------------ Initial abundances ------------
Y_CO = ab.Y_CO_equil(par.Y_C_tot, par.Y_O_tot, n, K_ra, K_rd)
X_CO = ph.A_CO*Y_CO
X_C_free = par.X_C_tot - X_CO
X_O_free = par.X_O_tot - X_CO
Y_C_free = X_C_free/ph.A_C
Y_O_free = X_O_free/ph.A_O
# Y_dust refers to the number fraction of carbon molecules in solid form
Y_solidC = 0.
n_solidC = Y_solidC*n
J_C = 0.
ncritical = 0.
# mass fractions are just A*Y, but dust mass fraction depends on grain size
X_solidC = ph.A_C*Y_solidC
# ---- CHECK: breaking up solid C mass frac into two components to see which one grows
X_dust_sml = 0.
X_dust_lrg = 0.
# ---------
n_C_free = Y_C_free * n
n_O_free = Y_O_free * n
n_CO = Y_CO * n
adap_flag = 0 # flag = 1 when dt is smaller than default


#------------ Initial dust parameters ------------
# size of grains:
# initialize a large array, and make a bigger one later if you need to
size = np.zeros([par.arraylen])
# similarly, save the number of grains at each time: 
Y_grains_C = np.zeros([par.arraylen])
# number fraction of grains at all times:
Y_grains_C_tot = np.zeros([par.arraylen])
neq=6.9e13*np.exp(-84428.2/T_cs)/ph.kB/T_cs
if neq==0: neq=par.floor   # just set it to a very small number if it comes out as 0
if neq!=neq: neq=par.floor # or if it comes out to be a NaN
# saturation:
sat_C = n_C_free/neq


# Initial time
t = par.tmin

#------------- OUTPUT FILES -----------------------#
# Create the run directory:
if not os.path.exists('runs/%s'%par.directory): os.makedirs('runs/%s'%par.directory)
#if os.path.exists("runs/%s/"%par.directory): print 'Directory exists, overwrite?'
# First open to overwrite previous - probably a better way to do this
fractionfile = open("runs/%s/fractions.txt"%par.directory, 'w')
densfile = open("runs/%s/densities.txt"%par.directory, 'w')
ratesfile = open("runs/%s/rates.txt"%par.directory, 'w')
thermofile = open("runs/%s/thermo.txt"%par.directory, 'w')
# Then open for appending
fractionfile = open("runs/%s/fractions.txt"%par.directory, 'a')
densfile = open("runs/%s/densities.txt"%par.directory, 'a')
ratesfile = open("runs/%s/rates.txt"%par.directory, 'a')
thermofile = open("runs/%s/thermo.txt"%par.directory, 'a')
# headers: 
fractionfile.write("# t     dt     X_CO      X_C_free     X_O_free    X_dust     X_D_small   X_D_large  int_flag    adap_flag     sat_C\n")
densfile.write("# t     n_CO      n_C_free     n_O_free     n_dust\n")
ratesfile.write("# K_ra      K_therm     K_nontherm   J_C \n")
thermofile.write("# T_cs     n\n")

# first line: 
fractionfile.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %i %.5e \n"%(t,par.dt_init,X_CO,X_C_free,X_O_free,X_solidC,X_dust_sml,X_dust_lrg,adap_flag,sat_C))
densfile.write("%.5f %.5f %.5f %.5f %.5f %.5f \n"%(t,n_CO,n_C_free,n_O_free,n_solidC,ncritical))
ratesfile.write("%.5f %.5e %.5e %.5e \n"%(K_ra, K_th, K_nth, J_C))
thermofile.write("%.5f %.5f\n"%(T_cs, n))


#----------- TIME EVOLUTION ----------------------# 
i = 1
# running to tmax - dt because interpolate function doesn't like t values outside of t_array
while t < par.tmax - par.dt_init:
	#---------------- First choose dt ------------------#
	dt = par.dt_init

	# -------------- Molecular Abundances ------------------
	#adap_flag, dt = choosedt(dt,n,neq,T_cs,Y_CO,Y_C_free,Y_O_free,K_ra,K_rd,n_C_free,size,dust_Y)

	# Evolve Y_CO
	dYCOdt = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd)		

	# ----------------------- DUST -----------------------

	# The equilibrium number fraction of carbon at this T
	neq=6.9e13*np.exp(-84428.2/T_cs)/ph.kB/T_cs
	if neq==0: neq=par.floor   # just set it to a very small number if it comes out as 0
	if neq!=neq: neq=par.floor # or if it comes out to be a NaN
	# Carbon saturation
	sat_C = n_C_free/neq 
	# Reset the rate of gaseous C depletion due to nucleation as 0
	dYdt_nucl_C = 0 
	J_C = 0.
	ncritical = 0
	# Compute the growth rate of previously formed grains [# of atoms / s]
	# depends on current size of each set of grains 
	dNdt_grow = np.zeros([par.arraylen])
	'''NO GROWTH FOR TEST'''
	#if sat_C > 1:
	#	dNdt_grow = ph.cshape*(size*ph.vC)**(2./3.) * np.sqrt(ph.kB*T_cs/(2*np.pi*ph.mC)) * n_C_free 
	# and update their sizes
	'''This also needs to move till after update time step'''
	size = size +dNdt_grow*dt

	# Indices of where the small/large grains are for summing in two bins later
	indx_sml = np.where(size < 1e6)
	indx_lrg = np.where(size > 1e6)

	# exception handling:
	# if the array gets filled up, break.
	try: 
		Y_grains_C[i] = 0.
	except IndexError:
		print 'Arrays are full - Exiting at t = ', t
		break

	if sat_C > 1: 
		# Compute nucleation rate [number of grains / cm^3 / s]
		J_C = (ph.cshape**3 * ph.vC**2 * ph.sigma/(18*np.pi**2*ph.mC))**(1./2.) * n_C_free**2 * \
			np.exp(-4*ph.cshape**3*ph.vC**2*ph.sigma**3/27./(ph.kB*T_cs)**3/(np.log(sat_C)**2)) 
		dmu = ph.kB*T_cs*np.log(sat_C)
		# Critical size [# of C atoms per grain]
		ncritical=np.maximum(8*ph.cshape**3 * ph.vC**2 * ph.sigma**3 /(27 * dmu**3),2)
		# Rate of depletion of gas C atoms due to nucleating new grains [number density of atoms / s]
		# AKA J*Ncrit*dt = number of carbon atoms that went into dust in this step
		dYdt_nucl_C = J_C*ncritical/n
		# Number of grains formed at each step: 
		Y_grains_C[i] = J_C*dt/n
		# initial size of newly formed grains [# of atoms]
		size[i]=ncritical


	# Add previous grains growth here!!

	# --------- Sum up some dust quantities for this time step ---------------
	# number fraction of present dust grains at current step
	''' this doesn't make sense unless it is multiplied by 
	 a Volume factor at each step (Vol is increasing) ''' 
	Y_grains_C_tot[i] = Y_grains_C.sum()

	# free C atoms that go into previous dust grain growth during this time step: 
	# Need to sum by number fraction
	'''check n factor here'''
	dYdt_growint_C = np.sum(dNdt_grow*Y_grains_C)



	# ------------- REFINE TIME STEP -------------------------
	delta_YC_CO = dYCOdt * dt
	delta_YC_nucl = dYdt_nucl_C*dt 
	delta_YC_grow = dYdt_growint_C*dt
	delta_YC_tot = (dYCOdt + dYdt_nucl_C + dYdt_growint_C)*dt

	while ((np.absolute(delta_YC_tot) > par.percent*Y_C_free + par.Y_min) \
		or np.absolute(delta_YC_tot) > par.percent*Y_O_free + par.Y_min):
		print 'timestep reduced ', Y_C_free, delta_YC_CO, delta_YC_nucl, delta_YC_grow
		adap_flag = 1
		dt = dt/2.
		delta_YC_tot = (dYCOdt + dYdt_nucl_C + dYdt_growint_C)*dt

	# --------------------------------------------------------

	t = t + dt


	# Y_solicC is the number fraction of carbon atoms that are locked in carbon clusters
	Y_solidC = Y_solidC + (dYdt_nucl_C+dYdt_growint_C)*dt
	X_solidC = ph.A_C*Y_solidC
	n_solidC = Y_solidC*n

	# Ok, not sure if this is the "correct" X because it doesn't correctly account for expansion,
	# but these X_D's are the total number of mass in either small grains or large grains,
	# we can see how these change over time
	X_dust_sml = ph.A_C*np.sum(Y_grains_C[indx_sml]*size[indx_sml])
	X_dust_lrg = ph.A_C*np.sum(Y_grains_C[indx_lrg]*size[indx_lrg])

	# --------------------------------------------------------


	# ----------------- Update the abundances ----------------
	# Y: number fraction
	# X: mass fraction
	# n: number density

	# CO 
	Y_CO = np.maximum(Y_CO +(dYCOdt)*dt,par.floor)
	X_CO = ph.A_CO*Y_CO
	n_CO = Y_CO * n
	# C goes into CO or nucleation or growth of carbon clusters
	Y_C_free = np.maximum(Y_C_free - delta_YC_tot,par.floor)
	Y_C_free = np.minimum(Y_C_free, par.Y_C_tot)
	X_C_free = ph.A_C*Y_C_free
	n_C_free = Y_C_free*n
	# O goes into CO
	Y_O_free = np.maximum(Y_O_free - (dYCOdt)*dt,par.floor)
	Y_O_free = np.minimum(Y_O_free, par.Y_O_tot)
	X_O_free = ph.A_O*Y_O_free	
	n_O_free = Y_O_free * n
	# -------------------------------------------------------
	if X_C_free != X_C_free: sys.exit('NaN!')
 
	# ---------- Append to text file ----------
	fractionfile.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %i %.5e \n"%(t,dt,X_CO,X_C_free,X_O_free,X_solidC,X_dust_sml,X_dust_lrg,adap_flag,sat_C))
	densfile.write("%.5f %.5f %.5f %.5f %.5f %.5f \n"%(t,n_CO,n_C_free,n_O_free,n_solidC,ncritical))
	ratesfile.write("%.5e %.5e %.5e %.5e \n"%(K_ra, K_th, K_nth, J_C))
	thermofile.write("%.5f %.5f \n"%(T_cs, n))

	# ------------ Expand the shell ------------
	# Evolve n, T, c_s, which here depend on R_cs, delta
	#R_cs, delta, n, T_cs, c_s = traj.shellexp(t, R_cs, delta, n, T_cs, c_s)
	#R_cs, n, T_cs = traj.pontefract(t,dt,R_cs,n,T_cs)
	n = n_interp(t)
	T_cs = T_interp(t)
	# ------------------------------------------

	# ------------ Evolve the rates ------------ 
	K_ra = re.formation(T_cs)
	K_th = re.thermal(T_cs)
	# Need to interpolate L_shock for non-thermal component
	# f_CO is the fraction of electrons that impact CO
	f_CO = Y_CO/(Y_CO+Y_C_free+Y_O_free+Y_solidC+par.Y_other)
	K_nth = L_interp(t)*par.epsilon*ph.mp*f_CO/par.M_cs/par.W_d
	K_rd = K_th + K_nth
	# ------------------------------------------

	#if X_C_free + X_O_free + X_CO +X_dust != 1.: 
	#	print 'ERROR: Sum of mass fractions is ne 1!'
	#	break

	if (int(t) % 10000100 < par.dt_init): print 't = ', t

	i = i+1

#---------------------------------------------------------#

# Save the 1d arrays for dust
dust_array = np.array([Y_grains_C, Y_grains_C_tot, size]).T
np.savetxt("runs/%s/dust.txt"%par.directory, dust_array, '%.5e', delimiter='   ',\
	header="# dust         alldust         size ")         


fractionfile.close()
ratesfile.close()
thermofile.close()
densfile.close()

