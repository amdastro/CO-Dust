import numpy as np
import pylab as plt
import physical as ph
import reaction as re
import parameters as par
import abundances as ab
import trajectory as traj
import os
import sys

# Create the run directory:
if not os.path.exists('runs/%s'%par.directory): os.makedirs('runs/%s'%par.directory)

#----------- Initial fluid parameters ------------
T_cs = par.T_cs_init
R_cs = par.R_cs_init
n = par.n_init
# For cold shell model: 
#c_s = np.sqrt(ph.kB * par.T_cs_init/(ph.mp))
#delta = par.delta_init

#------------ Initial reaction rates ------------
# In the future K_nonthermal will depend on time, but 
# for now it is constant.
K_nth = par.epsilon*par.L_gamma*ph.mp/par.M_cs/par.W_d
K_ra = re.formation(T_cs)
K_th = re.thermal(T_cs)
# K_rd is the sum of both destruction rates.
K_rd = K_th + K_nth

#------------ Initial abundances ------------
Y_CO = ab.Y_CO_equil(par.Y_C_tot, par.Y_O_tot, n, K_ra, K_rd)
X_CO = ph.A_CO*Y_CO
X_C_free = par.X_C_tot - X_CO
X_O_free = par.X_O_tot - X_CO
Y_C_free = X_C_free/ph.A_C
Y_O_free = X_O_free/ph.A_O
# Y_dust refers to the number fraction of carbon molecules in solid form
Y_dust = 0.
n_dust = Y_dust*n
J = 0.
ncritical = 0.
# mass fractions are just A*Y, but dust mass fraction depends on grain size
X_dust = ph.A_C*Y_dust
# ---- CHECK: breaking up solid C mass frac into two components to see which one grows
X_dust_sml = 0.
X_dust_lrg = 0.
# ---------
n_C_free = Y_C_free * n
n_O_free = Y_O_free * n
int_flag = 0 # flag = 1 when integrating
adap_flag = 0 # flag = 1 when dt is smaller than default

#------------ Initial dust parameters ------------
arraylen = 200000
# size of grains:
# initialize a large array, and make a bigger one later if you need to
size = np.zeros([arraylen])
# similarly, save the number of grains at each time: 
dust_Y = np.zeros([arraylen])
# number fraction of grains at all times:
alldust = np.zeros([arraylen])
# saturation:
sat = 0.


# Initial time
t = par.tmin

#------------- OUTPUT FILES -----------------------#
# First open to overwrite previous - probably a better way to do this
if os.path.exists("runs/%s/"%par.directory): print 'Directory exists, overwrite?'
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
fractionfile.write("# t     X_CO      X_C_free     X_O_free    X_dust     X_D_small   X_D_large  int_flag    adap_flag     sat\n")
densfile.write("# t     n_CO      n_C_free     n_O_free     n_dust\n")
ratesfile.write("# K_ra      K_therm     K_nontherm   J \n")
thermofile.write("# T_cs     n       R_cs   \n")
# first line: 
fractionfile.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %i %i %.5e \n"%(t,X_CO,X_C_free,X_O_free,X_dust,X_dust_sml,X_dust_lrg,int_flag,adap_flag,sat))
densfile.write("%.5f %.5f %.5f %.5f %.5f %.5f \n"%(t,n*Y_CO,n_C_free,n_O_free,n_dust,ncritical))
ratesfile.write("%.5f %.5e %.5e %.5e \n"%(K_ra, K_th, K_nth, J))
thermofile.write("%.5f %.5f %.5f \n"%(T_cs, n, R_cs))

#----------- TIME EVOLUTION ----------------------# 
i = 1
while t < par.tmax:
	#---------------- First choose dt ------------------#
	dt = par.dt_init
	# have function do all this:
	# adaptive_dt_flag, dt = timestep.choosedt(dt,n,T_cs, etc....)
	if t > par.t_integrate:
		int_flag = 1
		adap_flag = 0
		delta_YCO = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd) * dt
		''' 
		Should I change dt criteria for delta X_CO ???? 
		'''
		# Reduce dt if delta_YCO is too large
		while ((np.absolute(delta_YCO) > Y_C_free + par.CO_min) \
			or np.absolute(delta_YCO) > Y_O_free + par.CO_min):
			print 'reducing dt for Y_C_free ',Y_CO, delta_YCO, t 
			adap_flag = 1
			dt = dt/2.
			delta_YCO = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd) * dt
		# Also reduce dt so we don't form more than some percent of CO
		while (((Y_CO > par.CO_min) or (Y_CO + delta_YCO > par.CO_min)) \
			and (np.absolute(delta_YCO) > par.percentCO*Y_CO + par.CO_min)):
				print 'reducing dt for percent CO ',Y_CO, delta_YCO, t
				adap_flag = 1
				dt = dt/2.
				delta_YCO = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd) * dt

		# Also reduce dt for dust formation!
		# This is the prev step saturation and nucleation rate, but let's try this
		'''
		* * * *~**~~*~**~*~*~*~*~*~*~*~*~**~*~**~~*~*~*~*~*~*~* * * *
		 Need to create better time stepping criteria here 
		When dt is too large, create too much dust in first step and then mass is not conserved.
		* * * *~*~*~*~*~*~*~*~**~*~*~*~*****~~~~~*~**~*~*~*~*~* * * *
		'''
		sat = n_C_free/neq 
		if sat > 1:
			J = (ph.cshape**3 * ph.vC**2 * ph.sigma/(18*np.pi**2*ph.mC))**(1./2.) * n_C_free**2 * \
				np.exp(-4*ph.cshape**3*ph.vC**2*ph.sigma**3/27./(ph.kB*T_cs)**3/(np.log(sat)**2))
			dmu = ph.kB*T_cs*np.log(sat)
			ncritical=np.maximum(8*ph.cshape**3*ph.vC**2*ph.sigma**3/27/dmu**3,2)
			dYdt_nucl = J*ncritical/n
			delta_YC = dYdt_nucl*dt 
			'''Make a better timestepping criteria for dust formation!
			This doesn't include growth'''
			while (np.absolute(delta_YC) > 1e-2):#Y_C_free+par.CO_min):
				print 'reducing dt for dust ',dt
				dt = dt/2.
				delta_YC = dYdt_nucl*dt

		# Evolve Y_CO
		delta_YCO = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd) * dt
		Y_CO = np.maximum(delta_YCO + Y_CO,par.floor)
		delta_XCO = ph.A_CO*delta_YCO
		#X_CO = ph.A_CO*Y_CO
		# Subtract the number of atoms that went from free C and O to CO
		Y_C_free = np.maximum(Y_C_free - delta_YCO,par.floor)
		Y_O_free = np.maximum(Y_O_free - delta_YCO,par.floor)
		'''is this ok??? :'''
		Y_C_free = np.minimum(Y_C_free,par.Y_C_tot)
		Y_O_free = np.minimum(Y_O_free,par.Y_O_tot)
	else:
		# If not integrating, solve for the equilibrium density 
		int_flag = 0
		Y_CO = ab.Y_CO_equil(par.Y_C_tot, par.Y_O_tot, n, K_ra, K_rd)
		#X_CO = ph.A_CO*Y_CO
		Y_C_free = np.maximum(par.Y_C_tot - Y_CO,par.floor)
		Y_O_free = np.maximum(par.Y_O_tot - Y_CO,par.floor)

	#----------------------------------------------------#
	t = t + dt

	# update densities after CO evolution
	X_C_free = ph.A_C*Y_C_free
	X_O_free = ph.A_O*Y_O_free
	X_CO = ph.A_CO*Y_CO

	n_CO = Y_CO * n
	n_C_free = Y_C_free * n
	n_O_free = Y_O_free * n

	# dust: first compute the equilibrium number fraction of carbon at this T
	neq=6.9e13*np.exp(-84428.2/T_cs)/ph.kB/T_cs
	if neq==0: neq=par.floor   # just set it to a very small number if it comes out as 0
	if neq!=neq: neq=par.floor # or if it comes out to be a NaN
	# compute the saturation - free C over equil C
	sat = n_C_free/neq 
	# define the rate of gaseous C depletion due to nucleation as 0 (will change it later)
	dYdt_nucl = 0 
	J = 0.
	ncritical = 0
	# compute the growth rate of previously formed grains [# of atoms / s]
	# depends on current size of each set of grains 
	dNdt_grow = np.zeros([arraylen])
	if sat > 1:
		dNdt_grow = ph.cshape*(size*ph.vC)**(2./3.) * np.sqrt(ph.kB*T_cs/(2*np.pi*ph.mC)) * n_C_free 
	# and update their sizes
	size = size +dNdt_grow*dt

	# Indices of where the small/large grains are for summing in two bins later
	indx_sml = np.where(size < 1e6)
	indx_lrg = np.where(size > 1e6)

	# exception handling:
	# if the array gets filled up, break.
	try: 
		dust_Y[i] = 0.
	except IndexError:
		print 'Arrays are full - Exiting at t = ', t
		break

	if sat > 1: 
		# compute nucleation rate [number of grains / cm^3 / s]
		J = (ph.cshape**3 * ph.vC**2 * ph.sigma/(18*np.pi**2*ph.mC))**(1./2.) * n_C_free**2 * \
			np.exp(-4*ph.cshape**3*ph.vC**2*ph.sigma**3/27./(ph.kB*T_cs)**3/(np.log(sat)**2)) 
		dmu = ph.kB*T_cs*np.log(sat)
		# the critical size [# of C atoms per grain]
		ncritical=np.maximum(8*ph.cshape**3 * ph.vC**2 * ph.sigma**3 /(27 * dmu**3),2)
		# the rate of depletion of gas C atoms due to nucleating new grains [number density of atoms / s]
		# AKA J*Ncrit*dt = number of carbon atoms that went into dust in this step
		dYdt_nucl = J*ncritical/n
		# Number of grains formed at each step: 
		dust_Y[i] = J*dt/n
		# initial size of newly formed grains [# of atoms]
		size[i]=ncritical

	# number of present dust grains at current step
	''' this doesn't make sense unless it is multiplied by 
	 a Volume factor at each step (Vol is increasing) '''
	alldust[i] = dust_Y.sum()

	# free C atoms that go into previous dust grain growth during this time step: 
	# Need to sum by number fraction

	'''check n factor here'''
	dYdt_growint = np.sum(dNdt_grow*dust_Y)

	# Subtract from the free carbon mass fraction
	# Y_dust is the number fraction of carbon atoms that are in solid form
	Y_dust = Y_dust + (dYdt_nucl+dYdt_growint)*dt

	# Ok, not sure if this is the "correct" X because it doesn't correctly account for expansion,
	# but these X_D's are the total number of mass in either small grains or large grains,
	# we can see how these change over time

	X_dust_sml = ph.A_C*np.sum(dust_Y[indx_sml]*size[indx_sml])
	X_dust_lrg = ph.A_C*np.sum(dust_Y[indx_lrg]*size[indx_lrg])

	# Carbon dust, so 12 nucleons per atom:
	X_dust = ph.A_C*Y_dust
	n_dust = Y_dust*n


	# Save the mass distribution of grains at different times
	# save size array (number of atoms) and dust array (number of grains)
	# convert to mass dist in reading
	'''
	if (sat > 1 and t % 100000 < par.dt_init):
		massarray = np.array([dust_Y, size, dNdt_grow]).T
		print 'saving mass at t = ',t
		np.savetxt("runs/%s/massfile_t_%.0f.txt"%(par.directory,t), massarray, '%.5e', delimiter='   ',\
		header="# dust        size      growth rate")         
	'''

	# subtract the carbon from the free C gas
	Y_C_free = np.maximum(Y_C_free - (dYdt_nucl+dYdt_growint)*dt,par.floor)
	n_C_free = Y_C_free*n
	X_C_free = ph.A_C*Y_C_free

	# evolve the rates
	K_ra = re.formation(T_cs)
	K_th = re.thermal(T_cs)
	K_rd = K_th + K_nth
 
	# Append to text file
	fractionfile.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %i %i %.5e \n"%(t,X_CO,X_C_free,X_O_free,X_dust,X_dust_sml,X_dust_lrg,int_flag,adap_flag,sat))
	densfile.write("%.5f %.5f %.5f %.5f %.5f %.5f \n"%(t,n*Y_CO,n_C_free,n_O_free,n_dust,ncritical))
	ratesfile.write("%.5e %.5e %.5e %.5e \n"%(K_ra, K_th, K_nth, J))
	thermofile.write("%.5f %.5f %.5f \n"%(T_cs, n, R_cs))

	# Expand the shell
	# Evolve n, T, c_s, which here depend on R_cs, delta
	#R_cs, delta, n, T_cs, c_s = traj.shellexp(t, R_cs, delta, n, T_cs, c_s)
	R_cs, n, T_cs = traj.pontefract(t,dt,R_cs,n,T_cs)

	#if X_C_free + X_O_free + X_CO +X_dust != 1.: 
	#	print 'ERROR: Sum of mass fractions is ne 1!'
	#	break
	#print t
	i = i+1

#---------------------------------------------------------#

# Save the 1d arrays for dust
dust_array = np.array([dust_Y, alldust, size]).T
np.savetxt("runs/%s/dust.txt"%par.directory, dust_array, '%.5e', delimiter='   ',\
	header="# dust         alldust         size ")         


fractionfile.close()
ratesfile.close()
thermofile.close()
densfile.close()

