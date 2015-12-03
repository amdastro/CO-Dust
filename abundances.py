import physical as ph
import numpy as np


def n_CO_clayton(n_C, n_O, T) :
	'''Equilibrium density of CO from Clayton 2013 (to compare)'''
	nCO_clayton = n_C*n_O * (ph.h**2 / (2.*np.pi*ph.mu*ph.kB*T))**(1.5) * np.exp(ph.B_CO/(ph.kB*T))
	return nCO_clayton


def Y_CO_equil(Y_C, Y_O, n, Kra, Krd):
	'''Solution to quadratic eq for equilibrium number fraction of CO:
	Y_C and Y_O are the *total* number fractions.
	Although Y_C + Y_O = 1 here, I keep them separate since 
	this will not always be true. '''
	YCO_equil = 0.5*(Y_C + Y_O + Krd/Kra/n) - \
		0.5*((Y_C + Y_O + Krd/Kra/n)**2 - 4.*Y_C*Y_O)**0.5
	return YCO_equil


def dYCO_dt(Y_CO, Y_C, Y_O, n, Kra, Krd):
	'''Define function for the change in Y_CO
	this is a first order ODE
	Y_C and Y_O are the *free* number fractions here. '''
	dY_CO_dt = n*Kra*Y_C*Y_O - Krd*Y_CO
	return dY_CO_dt