import numpy as np
import pylab as plt
import physical as ph
import reaction as re
import parameters as par
import abundances as ab
import trajectory as traj
from scipy.interpolate import interp1d
import os
import sys

print 'RUNNING ', par.directory

#----------- Initial fluid parameters ------------
T_cs = par.T_cs_init
n = par.n_init
# Interpolated arrays from thermo data table
n_interp = interp1d(par.t_array,par.n_array,kind='linear')
T_interp = interp1d(par.t_array,par.T_array,kind='linear')
L_interp = interp1d(par.t_array,par.L_array,kind='linear')

#------------ Initial reaction rates ------------
K_CO_ra  = re.COformation(T_cs)
K_CO_th  = re.COthermal(T_cs)
K_SiO_ra = re.SiOformation(T_cs)
K_SiO_th = re.SiOthermal(T_cs)
# Start with the equilibrium mass fractions assuming
# there is no nonthermal destruction initially.
# Because we need the abundance to get the destruciton rate. 
K_CO_nth  = 0.
K_SiO_nth = 0.
K_CO_rd   = K_CO_th + K_CO_nth
K_SiO_rd  = K_SiO_th + K_SiO_nth


#------------ Initial abundances ------------
Y_CO  = ab.Y_CO_equil(par.Y_C_tot, par.Y_O_tot, n, K_CO_ra, K_CO_rd)
Y_SiO = ab.Y_SiO_equil(par.Y_Si_tot, par.Y_O_tot, Y_CO, n, K_SiO_ra, K_SiO_rd)
X_CO  = ph.A_CO*Y_CO
X_SiO = ph.A_SiO*Y_SiO
X_C   = par.X_C_tot - X_CO
Y_C   = X_C/ph.A_C
X_O   = par.X_O_tot - X_CO - X_SiO
Y_O   = X_O/ph.A_O
X_Si  = par.X_Si_tot - X_SiO
Y_Si  = X_Si / ph.A_Si
X_Mg  = par.X_Mg_tot
Y_Mg  = X_Mg / ph.A_Mg

n_C   = Y_C * n
n_O   = Y_O * n
n_CO  = Y_CO * n
n_SiO = Y_SiO * n
n_Si  = Y_Si * n
n_Mg  = Y_Mg * n


# Solid mass fractions
# Carbon atoms locked in grains
Y_solidC = 0.
n_solidC = Y_solidC * n
X_solidC = ph.A_C * Y_solidC
# Mg2SiO4 molecules locked in grains (won't sum to unity here!!)
Y_Mg2SiO4 = 0.
n_Mg2SiO4 = Y_Mg2SiO4 * n
# NEED to distinguish between Y_Mg2SiO4, which is the number fraction of existing Mg2SiO4 molcules
# Y_Sidust which is the number fraction of either Mg, (Si + O), or O that is in dust
X_Mg2SiO4 = ph.A_Mg2SiO4 * Y_Mg2SiO4
# Can make a "fake" number fraction that just adds the number of Si, O or Mg atoms that are in Si dust
# Be careful summing the mass fraction of this
# Then check to make sure mass is conserved. 

# flag = 1 when dt is smaller than default
adap_flag = 0 

#------------ Initial dust parameters ------------
# Carbon grains
J_Cg = 0.
ncritical_Cg = 0.
# Initialize large arrays, and make bigger ones later if you need to
# size of grains:
sizes_Cg = np.zeros([par.arraylen])
# save the number of grains at each time: 
Y_grains_Cg = np.zeros([par.arraylen])
# number fraction of grains at all times:
Y_grains_Cg_tot = np.zeros([par.arraylen])
# Initial saturation: 
neq=6.9e13*np.exp(-84428.2/T_cs)/ph.kB/T_cs
if neq==0: neq=par.floor   
if neq!=neq: neq=par.floor 
sat_Cg = n_C/neq

# Silicate grains
J_Sig = 0.
ncritical_Sig = 0.
# size of grains:
sizes_Sig = np.zeros([par.arraylen])
# save the number of grains at each time: 
Y_grains_Sig = np.zeros([par.arraylen])
# number fraction of grains at all times:
Y_grains_Sig_tot = np.zeros([par.arraylen])
# Initial saturation: 
# ln S = -G / (kT) + sum (nu_i * ln(p_i))
# partial pressures: Volume element?
'''CORREST PARTIAL PRESSURES'''
p_Mg = n_Mg*ph.kB*T_cs
p_SiO = n_SiO * ph.kB*T_cs
p_O = n_O * ph.kB*T_cs
lnsat_Sig = ph.GibbsA/T_cs - ph.GibbsB + (2.*np.log(p_Mg) + np.log(p_SiO) + 3.*np.log(p_O))


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
fractionfile.write("# t     dt     X_CO      X_C     X_O    X_solidC   \
	X_Mg      X_Si     X_SiO       X_Mg2SiO4    adap_flag     sat_Cg    lnsat_Sig\n")
densfile.write("# t     n_CO      n_C      n_O       n_solidC    n_Mg      \
	n_Si      n_SiO    n_Mg2SiO4  ncritical_Cg     ncritical_Sig\n")
ratesfile.write("# K_CO_ra      K_CO_therm     K_CO_nontherm   K_SiO_ra     \
	K_SiO_therm      K_SiO_nontherm     J_Cg     J_Sg\n")
thermofile.write("# T_cs     n\n")

# first line: 
fractionfile.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %i %.5f %.5e \n"%(t,par.dt_init,\
	X_CO,X_C,X_O,X_solidC,X_Mg,X_Si,X_SiO,X_Mg2SiO4,adap_flag,sat_Cg,lnsat_Sig))
densfile.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f \n"%(t,\
	n_CO,n_C,n_O,n_solidC,n_Mg,n_Si,n_SiO,n_Mg2SiO4,ncritical_Cg,ncritical_Sig))
ratesfile.write("%.5e %.5e  %.5e %.5e %.5e %.5e %.5e %.5e \n"%(K_CO_ra, \
	K_CO_th, K_CO_nth, K_SiO_ra, K_SiO_th, K_SiO_nth,J_Cg,J_Sig))
thermofile.write("%.5f %.5f\n"%(T_cs, n))


#----------- TIME EVOLUTION ----------------------# 
i = 1
# running to tmax - dt because interpolate function doesn't like t values outside of t_array
while t < par.tmax - par.dt_init:
	#---------------- First choose dt ------------------#
	dt = par.dt_init

	# -------------- Molecular Abundances ------------------
	# Evolve Y_CO
	dYCOdt = ab.dYCO_dt(Y_CO, Y_C, Y_O, n, K_CO_ra, K_CO_rd)		
	dYSiOdt = ab.dYSiO_dt(Y_SiO, Y_Si, Y_O, n, K_SiO_ra, K_SiO_rd)
	# ----------------------- DUST -----------------------

	# The equilibrium number fraction of carbon at this T
	neq=6.9e13*np.exp(-84428.2/T_cs)/ph.kB/T_cs
	if neq==0: neq=par.floor   # just set it to a very small number if it comes out as 0
	if neq!=neq: neq=par.floor # or if it comes out to be a NaN
	# Carbon saturation
	sat_Cg = n_C/neq 
	# Silicate saturation: 
	# ln S = -G / (kT) + sum (nu_i * p_i)
	# partial pressures: Volume element?
	p_Mg = n_Mg*ph.kB*T_cs
	p_SiO = n_SiO * ph.kB*T_cs
	p_O = n_O * ph.kB*T_cs
	lnsat_Sig = ph.GibbsA/T_cs - ph.GibbsB + (2.*np.log(p_Mg) + np.log(p_SiO) + 3.*np.log(p_O))


	# Reset the rate of gaseous C depletion due to nucleation
	dYdt_nucl_Cg = 0 
	J_Cg = 0.
	ncritical_Cg = 0
	# Reset the rate of gaseous Si,Mg,O depletion due to nucleation
	dYdt_nucl_Sig = 0 
	J_Sig = 0.
	ncritical_Sig = 0
	# Reset the growth rate of previously formed grains [# of atoms / s]
	dNdt_grow_Cg  = np.zeros([par.arraylen])
	dNdt_grow_Sig = np.zeros([par.arraylen])

	# exception handling:
	# if the array gets filled up, break.
	try: 
		Y_grains_Cg[i] = 0.
	except IndexError:
		print 'Arrays are full - Exiting at t = ', t
		break

	'''NO DUST'''
	sat_Cg = 0
	if sat_Cg > 1: 
		J_Cg = (ph.cshape**3 * ph.vC**2 * ph.sigmaC/(18*np.pi**2*ph.mC))**(1./2.) * n_C**2 * \
		  np.exp(-4*ph.cshape**3*ph.vC**2*ph.sigmaC**3/27./(ph.kB*T_cs)**3/(np.log(sat_Cg)**2)) 
		#'''NO GROWTH FOR TEST'''
		# Grow the previous grains
		# dNdt_grow_Cg = ph.cshape*(sizes_Cg*ph.vC)**(2./3.) * np.sqrt(ph.kB*T_cs/(2*np.pi*ph.mC)) * n_C 
		# Compute nucleation rate [number of grains / cm^3 / s]  
		dmu_Cg = ph.kB*T_cs*np.log(sat_Cg)
		# Critical size [# of C atoms per grain]
		ncritical_Cg = np.maximum(8*ph.cshape**3 * ph.vC**2 * ph.sigmaC**3 /(27 * dmu_Cg**3),2)
		# Rate of depletion of gas C atoms due to nucleating new grains [number density of atoms / s]
		# AKA J*Ncrit*dt = number of carbon atoms that went into dust in this step
		dYdt_nucl_Cg = J_Cg*ncritical_Cg/n
		# Number of grains formed at each step: 
		Y_grains_Cg[i] = J_Cg*dt/n
		# initial size of newly formed grains [# of atoms]
		sizes_Cg[i]=ncritical_Cg

	'''NO DUST'''
	lnsat_Sig = 0
	if lnsat_Sig > 0: 
		J_Sig = (ph.cshape**3 * ph.vMg2SiO4**2 * ph.sigmaMg2SiO4/(18*np.pi**2*ph.mMg2SiO4))**(1./2.) * n_C**2 * \
			np.exp(-4*ph.cshape**3*ph.vMg2SiO4**2*ph.sigmaMg2SiO4**3/27./(ph.kB*T_cs)**3/lnsat_Sig**2 )
		#'''NO GROWTH FOR TEST'''
		# Grow the previous grains
		# dNdt_grow_Sig = ph.cshape*(sizes_Sig*ph.vC)**(2./3.) * np.sqrt(ph.kB*T_cs/(2*np.pi*ph.mC)) * n_C 
		# Compute nucleation rate [number of grains / cm^3 / s]
		dmu_Sig = ph.kB*T_cs*lnsat_Sig
		# Critical size [# of Mg2SiO4 molecules per grain]
		ncritical_Sig = np.maximum(8*ph.cshape**3 * ph.vMg2SiO4**2 * ph.sigmaMg2SiO4**3 /(27 * dmu_Sig**3),2)
		# Rate of depletion of gas into Mg2SiO4 due to nucleating new grains [number density of molecules / s]
		# AKA J*Ncrit*dt = number of Mg2SiO4 molecules that went into dust in this step
		dYdt_nucl_Sig = J_Sig*ncritical_Sig/n
		# Number of grains formed at each step: 
		Y_grains_Sig[i] = J_Sig*dt/n
		# initial size of newly formed grains [# of atoms]
		sizes_Sig[i]=ncritical_Sig


	# --------- Sum up some dust quantities for this time step ---------------
	# number fraction of present dust grains at current step
	''' this doesn't make sense unless it is multiplied by 
	 a Volume factor at each step (Vol is increasing) ''' 
	Y_grains_Cg_tot[i] = Y_grains_Cg.sum()
	Y_grains_Sig_tot[i] = Y_grains_Sig.sum()

	# free C atoms that go into previous dust grain growth during this time step: 
	# Need to sum by number fraction
	'''check n factor here'''
	dYdt_growint_Cg = np.sum(dNdt_grow_Cg*Y_grains_Cg)
	dYdt_growint_Sig = np.sum(dNdt_grow_Sig*Y_grains_Sig)
	# Abundance changes from silicate grains
	dYMgdt_nucl  = 2.*dYdt_nucl_Sig
	dYSiOdt_nucl = dYdt_nucl_Sig
	dYOdt_nucl   = 3.*dYdt_nucl_Sig
	#dYMgdt_grow  = 2.*dNdt_grow_Sig / n 
	#dYSiOdt_grow = dNdt_grow_Sig / n
	#dYOdt_grow   = 3.*dNdt_grow_Sig / n
	dYMgdt_grow  = 2.*dYdt_growint_Sig
	dYSiOdt_grow = dYdt_growint_Sig
	dYOdt_grow   = 3.*dYdt_growint_Sig

	# ------------- REFINE TIME STEP -------------------------
	# Sum up changes from molecular reactions, nucleation, and growth
	delta_YC_tot  = (-dYCOdt - dYdt_nucl_Cg - dYdt_growint_Cg)*dt
	delta_YO_tot  = (-(dYCOdt + dYSiOdt) - dYOdt_nucl - dYdt_growint_Sig) * dt
	delta_YMg_tot = (-dYMgdt_nucl - dYMgdt_grow) * dt
	delta_YSi_tot = (-dYSiOdt) * dt
	delta_YSiO_tot= (dYSiOdt - dYSiOdt_nucl - dYSiOdt_grow) * dt

	# right now the dt criteria is only set for C and O

	while ((np.absolute(delta_YC_tot) > par.percent*Y_C + par.Y_min) \
		or (np.absolute(delta_YO_tot) > par.percent*Y_O + par.Y_min)):
		print 'timestep reduced '
		adap_flag = 1
		dt = dt/2.
		delta_YC_tot  = (-dYCOdt - dYdt_nucl_Cg - dYdt_growint_Cg)*dt
		delta_YO_tot  = (-(dYCOdt + dYSiOdt) - dYOdt_nucl - dYdt_growint_Sig) * dt
	
	# --------------------------------------------------------
	# Update grain sizes
	sizes_Cg = sizes_Cg +dNdt_grow_Cg*dt
	sizes_Sig = sizes_Sig +dNdt_grow_Sig*dt
	# Update time
	t = t + dt


	# Y_solicC is the number fraction of carbon atoms that are locked in carbon clusters
	Y_solidC = Y_solidC + (dYdt_nucl_Cg+dYdt_growint_Cg)*dt
	X_solidC = ph.A_C*Y_solidC
	n_solidC = Y_solidC*n

	#Note: this mass fraction will not add to 1
	Y_Mg2SiO4 = Y_Mg2SiO4 + (dYdt_nucl_Sig+dYdt_growint_Sig)*dt
	n_Mg2SiO4 = Y_Mg2SiO4 * n
	X_Mg2SiO4 = ph.A_Mg2SiO4 * Y_Mg2SiO4
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
	Y_C = np.maximum(Y_C + delta_YC_tot,par.floor)
	Y_C = np.minimum(Y_C, par.Y_C_tot)
	X_C = ph.A_C*Y_C
	n_C = Y_C*n
	# O goes into CO and SiO and Mg2SiO4
	Y_O = np.maximum(Y_O + delta_YO_tot,par.floor)
	Y_O = np.minimum(Y_O, par.Y_O_tot)
	X_O = ph.A_O*Y_O	
	n_O = Y_O * n
	# Mg goes into Mg2SiO4
	Y_Mg = np.maximum(Y_Mg + delta_YMg_tot, par.floor)
	X_Mg = ph.A_Mg*Y_Mg
	n_Mg = Y_Mg * n
	# Si goes into SiO
	Y_Si = np.maximum(Y_Si + delta_YSi_tot, par.floor)
	Y_Si = np.minimum(Y_Si, par.Y_Si_tot)
	X_Si = ph.A_Si*Y_Si
	n_Si = Y_Si * n
	# SiO goes into Mg2SiO4
	Y_SiO = np.maximum(Y_SiO + delta_YSiO_tot, par.floor)
	X_SiO = ph.A_SiO*Y_SiO
	n_SiO = Y_SiO * n

	# -------------------------------------------------------
	if X_C != X_C: sys.exit('NaN!')
 
	# ---------- Append to text file ----------
	fractionfile.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %i %.5f %.5e \n"%(t,par.dt_init,\
		X_CO,X_C,X_O,X_solidC,X_Mg,X_Si,X_SiO,X_Mg2SiO4,adap_flag,sat_Cg,lnsat_Sig))
	densfile.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f \n"%(t,\
		n_CO,n_C,n_O,n_solidC,n_Mg,n_Si,n_SiO,n_Mg2SiO4,ncritical_Cg,ncritical_Sig))
	ratesfile.write("%.5e %.5e  %.5e %.5e %.5e %.5e %.5e %.5e \n"%(K_CO_ra,\
	 K_CO_th, K_CO_nth, K_SiO_ra, K_SiO_th, K_SiO_nth,J_Cg,J_Sig))
	thermofile.write("%.5f %.5f\n"%(T_cs, n))

	# ------------ Expand the shell ------------
	# Evolve n, T, c_s, which here depend on R_cs, delta
	#R_cs, delta, n, T_cs, c_s = traj.shellexp(t, R_cs, delta, n, T_cs, c_s)
	#R_cs, n, T_cs = traj.pontefract(t,dt,R_cs,n,T_cs)
	n = n_interp(t)
	T_cs = T_interp(t)
	# ------------------------------------------

	# ------------ Evolve the rates ------------ 
	K_CO_ra  = re.COformation(T_cs)
	K_CO_th  = re.COthermal(T_cs)
	K_SiO_ra = re.SiOformation(T_cs)
	K_SiO_th = re.SiOthermal(T_cs)
	# Need to interpolate L_shock for non-thermal component
	# f_A is the fraction of electrons that impact element A
	f_CO = Y_CO/(Y_CO+Y_C+Y_O+Y_Si+Y_Mg+Y_SiO+par.Y_other)
	f_SiO = Y_SiO/(Y_CO+Y_C+Y_O+Y_Si+Y_Mg+Y_CO+par.Y_other)
	K_CO_nth = L_interp(t)*par.epsilon*ph.mp*f_CO/par.M_cs/par.W_d_CO
	K_SiO_nth = L_interp(t)*par.epsilon*ph.mp*f_SiO/par.M_cs/par.W_d_SiO
	K_CO_rd = K_CO_th + K_CO_nth
	K_SiO_rd = K_SiO_th + K_SiO_nth
	# ------------------------------------------

	if (int(t) % 10000100 < par.dt_init): print 't = ', t
	i = i+1

#---------------------------------------------------------#

# Save the 1d arrays for dust
Cdust_array = np.array([Y_grains_Cg, Y_grains_Cg_tot, sizes_Cg]).T
np.savetxt("runs/%s/Cdust.txt"%par.directory, Cdust_array, '%.5e', delimiter='   ',\
	header="# dust         alldust         size ")         

Sidust_array = np.array([Y_grains_Sig, Y_grains_Sig_tot, sizes_Sig]).T
np.savetxt("runs/%s/Sidust.txt"%par.directory, Sidust_array, '%.5e', delimiter='   ',\
	header="# dust         alldust         size ")   

fractionfile.close()
ratesfile.close()
thermofile.close()
densfile.close()

