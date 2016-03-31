import physical as ph
import numpy as np



def Y_CO_equil(Y_C, Y_O, n, Kra, Krd):
	'''Solution to quadratic eq for equilibrium number fraction of CO:
	Y_C and Y_O are the *total* number fractions.'''
	YCO_equil = 0.5*(Y_C + Y_O + Krd/Kra/n) - \
		0.5*((Y_C + Y_O + Krd/Kra/n)**2 - 4.*Y_C*Y_O)**0.5
	return YCO_equil


def dYCO_dt(Y_CO, Y_C, Y_O, n, Kra, Krd):
	'''Define function for the change in Y_CO
	this is a first order ODE
	Y_C and Y_O are the *free* number fractions here. '''
	dY_CO_dt = n*Kra*Y_C*Y_O - Krd*Y_CO
	return dY_CO_dt


def Y_SiO_equil(Y_Si, Y_O, Y_CO, n, Kra, Krd):
	'''Solution to quadratic eq for equilibrium number fraction of SiO:
	Y_Si and Y_O are the *total* number fractions.'''
	#YSiO_equil = 0.5*(Y_Si + Y_O + Krd/Kra/n) - \
	#	0.5*((Y_Si + Y_O + Krd/Kra/n)**2 - 4.*Y_Si*Y_O)**0.5
	YSiO_equil = 0.5*(Krd/(n*Kra) - Y_CO + Y_O + Y_Si) - \
		0.5*(4*Y_Si*(Y_CO - Y_O) +(Krd + n*Kra*(-Y_CO + Y_O + Y_Si))**2/Kra**2/n**2)**0.5
	return YSiO_equil


def dYSiO_dt(Y_SiO, Y_Si, Y_O, n, Kra, Krd):
	'''Define function for the change in Y_CO
	this is a first order ODE
	Y_C and Y_O are the *free* number fractions here. 
	Rates are different than for CO'''
	dYSiO_dt = n*Kra*Y_Si*Y_O - Krd*Y_SiO
	return dYSiO_dt