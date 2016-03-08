import numpy as np
import abundances as ab
import physical as ph
import parameters as par

def choosedt(dt,n,neq,T_cs,Y_CO,Y_C_free,Y_O_free,K_ra,K_rd,n_C_free,size,dust_Y):
	adap_flag = 0
	# Change in Y_C_free from CO formation:
	delta_YC_CO = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd) * dt
	# Change in Y_C_free from carbon dust nucleation (and growth)
	dYdt_nucl = 0.
	dNdt_grow = np.zeros([par.arraylen])
	sat = n_C_free/neq 
	if sat > 1:
		J = (ph.cshape**3 * ph.vC**2 * ph.sigma/(18*np.pi**2*ph.mC))**(1./2.) * n_C_free**2 * \
			np.exp(-4*ph.cshape**3*ph.vC**2*ph.sigma**3/27./(ph.kB*T_cs)**3/(np.log(sat)**2))
		dmu = ph.kB*T_cs*np.log(sat)
		ncritical=np.maximum(8*ph.cshape**3*ph.vC**2*ph.sigma**3/27/dmu**3,2)
		dYdt_nucl = J*ncritical/n	
		# Calculate the possible growth
		dNdt_grow = ph.cshape*(size*ph.vC)**(2./3.) * np.sqrt(ph.kB*T_cs/(2*np.pi*ph.mC)) * n_C_free
		delta_YC_grow = np.sum(dNdt_grow*dust_Y)*dt # dust_Y is a 1d array storing previous dYdt_nucl[i], read into function 

	delta_YC_nucl = dYdt_nucl*dt 
	delta_YC_grow = np.sum(dNdt_grow*dust_Y)*dt

	# Maybe this will be too stringent, but let's try
	# if the sum (not absolute value?) of all changes in Y_C_free is greater than 
	# say, 30% of Y_C_free (+ a minimum value) or Y_O_free, for future abundances
	# then reduce dt:
	delta_YC_tot = delta_YC_CO + delta_YC_nucl + delta_YC_grow
	while ((np.absolute(delta_YC_tot) > par.percent*Y_C_free + par.Y_min) \
		or np.absolute(delta_YC_tot) > par.percent*Y_O_free + par.Y_min):
		print 'timestep reduced ', Y_C_free, delta_YC_CO, delta_YC_nucl, delta_YC_grow
		adap_flag = 1
		dt = dt/2.
		delta_YC_CO = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd) * dt
		delta_YC_nucl = dYdt_nucl*dt
		delta_YC_grow = np.sum(dNdt_grow*dust_Y)*dt
		delta_YC_tot = delta_YC_CO + delta_YC_nucl + delta_YC_grow

	return adap_flag, dt




