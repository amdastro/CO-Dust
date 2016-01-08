import numpy as np
import pylab as plt
import physical as ph
import reaction as re
import parameters as par
import abundances as ab
import os
import sys

# Create the run directory:
if not os.path.exists('runs/%s'%par.directory): os.makedirs('runs/%s'%par.directory)

#----------- Initial fluid parameters ------------
T_cs = par.T_cs_init
c_s = np.sqrt(ph.kB * par.T_cs_init/(ph.mp))
R_cs = par.R_cs_init
delta = par.delta_init
n = par.n_init

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
# mass fractions are just A*Y, but dust mass fraction depends on grain size
X_dust = ph.A_C*Y_dust
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
dust = np.zeros([arraylen])
# number of grains at all times:
alldust = np.zeros([arraylen])
# saturation:
sat = 0.
# all carbon in solid form: 
allcarbon = np.zeros([arraylen])

# Initial time
t = par.tmin

#------------- OUTPUT FILES -----------------------#
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
fractionfile.write("# t     X_CO      X_C_free     X_O_free    X_dust    int_flag    adap_flag     sat\n")
densfile.write("# t     n_CO      n_C_free     n_O_free     n_dust\n")
ratesfile.write("# K_ra      K_therm     K_nontherm \n")
thermofile.write("# T_cs     n      delta      R_cs       c_s  \n")
# first line: 
fractionfile.write("%.5f %.5f %.5f %.5f %.5f %i %i %.5e \n"%(t,X_CO,X_C_free,X_O_free,X_dust,int_flag,adap_flag,sat))
densfile.write("%.5f %.5f %.5f %.5f %.5f \n"%(t,n*Y_CO,n_C_free,n_O_free,n_dust))
ratesfile.write("%.5f %.5e %.5e \n"%(K_ra, K_th, K_nth))
thermofile.write("%.5f %.5f %.5f %.2f %.5f \n"%(T_cs, n, delta, R_cs, c_s))

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
		sat = n_C_free/neq 
		if sat > 1:
			J = (ph.cshape**3 * ph.vC**2 * ph.sigma/(18*np.pi**2*ph.mC))**(1./2.) * n_C_free**2 * \
				np.exp(-4*ph.cshape**3*ph.vC**2*ph.sigma**3/27./(ph.kB*T_cs)**3/(np.log(sat)**2))
			dmu = ph.kB*T_cs*np.log(sat)
			ncritical=np.maximum(8*ph.cshape**3*ph.vC**2*ph.sigma**3/27/dmu**3,2)
			dCdt_nucl = J*ncritical
			delta_nC = dCdt_nucl*dt
			'''Make a better timestepping criteria for dust formation!'''
			while (np.absolute(delta_nC) > n_C_free+par.CO_min):
				print 'reducing dt for dust ',dt
				dt = dt/2.
				delta_dust = dYdust_dt_nucl*dt

		# Evolve Y_CO
		delta_YCO = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd) * dt
		Y_CO = np.maximum(delta_YCO + Y_CO,1e-100)
		delta_XCO = ph.A_CO*delta_YCO
		#X_CO = ph.A_CO*Y_CO
		# Subtract the number of atoms that went from free C and O to CO
		Y_C_free = np.maximum(Y_C_free - delta_YCO,1e-100)
		Y_O_free = np.maximum(Y_O_free - delta_YCO,1e-100)
		'''is this ok??? :'''
		Y_C_free = np.minimum(Y_C_free,par.Y_C_tot)
		Y_O_free = np.minimum(Y_O_free,par.Y_O_tot)
	else:
		# If not integrating, solve for the equilibrium density 
		int_flag = 0
		Y_CO = ab.Y_CO_equil(par.Y_C_tot, par.Y_O_tot, n, K_ra, K_rd)
		#X_CO = ph.A_CO*Y_CO
		Y_C_free = np.maximum(par.Y_C_tot - Y_CO,1e-100)
		Y_O_free = np.maximum(par.Y_O_tot - Y_CO,1e-100)

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
	if neq==0: neq=1e-280   # just set it to a very small number if it comes out as 0
	if neq!=neq: neq=1e-280 # or if it comes out to be a NaN
	# compute the saturation - free C over equil C
	sat = n_C_free/neq 
	# define the rate of gaseous C depletion due to nucleation as 0 (will change it later)
	dCdt_nucl = 0 
	# compute the growth rate of previously formed grains [# of atoms / s]
	# depends on current size of each set of grains
	dNdt_grow = ph.cshape*(size[:i-1]*ph.vC)**(2./3.) * np.sqrt(ph.kB*T_cs/(2*np.pi*ph.mC)) * n_C_free 
	# and update their sizes
	size[:i-1] = size[:i-1]+dNdt_grow[:i-1]*dt

	# exception handling:
	# if the array gets filled up, break
	try: 
		dust[i] = 0.
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
		dCdt_nucl = J*ncritical
		# Number of grains formed at each step: 
		dust[i] = J*dt
		# initial size of newly formed grains [# of atoms]
		size[i]=ncritical

	# number density of present dust grains at current step
	alldust[i] = dust[:i].sum()
	# total present solid carbon atoms at current step
	''' Not a true number density because it doesn't take into account the expansion! '''
	allcarbon[i] = np.sum(dust[:i]*size[:i])

	# free C atoms that go into previous dust grain growth during this time step: 
	dCdt_growint = np.sum(dNdt_grow[:i-1]*dust[:i-1])

	# Subtract from the free carbon mass fraction
	# Y_dust is the number fraction of carbon atoms that are in solid form
	Y_dust = np.minimum(Y_dust + (dCdt_nucl+dCdt_growint)/n*dt,par.Y_C_tot)
	# Carbon dust, so 12 nucleons per atom:
	X_dust = ph.A_C*Y_dust
	n_dust = Y_dust*n

	# subtract the carbon from the free C gas
	n_C_free = np.maximum(n_C_free - (dCdt_nucl+dCdt_growint)*dt,0)
	Y_C_free = n_C_free/n
	X_C_free = ph.A_C*Y_C_free

	# evolve the rates
	K_ra = re.formation(T_cs)
	K_th = re.thermal(T_cs)
	K_rd = K_th + K_nth
 
	# Append to text file
	fractionfile.write("%.5f %.5f %.5f %.5f %.5f %i %i %.5e \n"%(t,X_CO,X_C_free,X_O_free,X_dust,int_flag,adap_flag,sat))
	densfile.write("%.5f %.5f %.5f %.5f %.5f \n"%(t,n*Y_CO,n_C_free,n_O_free,n_dust))
	ratesfile.write("%.5e %.5e %.5e \n"%(K_ra, K_th, K_nth))
	thermofile.write("%.5f %.5f %.5f %.2f %.5f \n"%(T_cs, n, delta, R_cs, c_s))

	# Expand the shell
	R_cs = par.R_cs_init + par.v_ej*t
	delta = par.delta_init + c_s*t
	n = par.M_cs/(4.*np.pi*ph.mp*R_cs**2 * delta)

	# Temperature and sound speed decrease
	T_cs = par.T_cs_init * (n / par.n_init)**(par.gamma-1.)
	c_s = np.sqrt(ph.kB*T_cs/(ph.mp))

	#if X_C_free + X_O_free + X_CO +X_dust != 1.: 
	#	print 'ERROR: Sum of mass fractions is ne 1!'
	#	break
	#print t
	i = i+1
#---------------------------------------------------------#

# Save the 1d arrays for dust
dust_array = np.array([dust, alldust, size, allcarbon]).T
np.savetxt("runs/%s/dust.txt"%par.directory, dust_array, '%.5e', delimiter='   ',\
	header="# dust         alldust         size         allcarbon")         


fractionfile.close()
ratesfile.close()
thermofile.close()
densfile.close()

