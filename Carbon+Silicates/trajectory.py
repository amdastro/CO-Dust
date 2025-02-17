import numpy as np
import parameters as par
import physical as ph



def shellexp(t,R_cs,delta,n,T_cs,c_s):
	R_cs = par.R_cs_init + par.v_ej*t
	delta = par.delta_init + c_s*t
	n = par.M_cs/(4.*np.pi*ph.mp*R_cs**2 * delta)
	T_cs = par.T_cs_init * (n / par.n_init)**(par.gamma-1.)
	c_s = np.sqrt(ph.kB*T_cs/(ph.mp))
	return R_cs,delta,n,T_cs,c_s


def pontefract(t,dt,R,R0,v_ej,n,n0,T,T0):
	R = R + v_ej*dt
	n = n0 * (R/R0)**-3
	T = T0 * (R/R0)**-0.5
	return R,n,T


#def shocklum(t):
	'''Read the next temperature and density from table
	Might not do this in function but in main code'''
#	index = np.where(par.t_array > t & par.t_array <= t+dt)
#	print 't, index, ti-1, ti+1 = ', t, index
#	n = n_array[index[0]]
#	T_cs = T_array[index[0]]
#	return n, T_cs