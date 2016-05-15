import physical as ph
import numpy as np



def Y_CO_equil(Y_C, Y_O, Y_SiO, n, Kra, Krd):
	'''Solution to quadratic eq for equilibrium number fraction of CO:
	Y_C and Y_O and Y_H are the *total* number fractions.'''
	#YCO_equil = 0.5*(Y_C + Y_O + Krd/Kra/n) - \
	#	0.5*((Y_C + Y_O + Krd/Kra/n)**2 - 4.*Y_C*Y_O)**0.5
	# Now adding SiO and CO+H reaction
	#YCO_equil = 0.5*(Y_O + Y_C - Y_SiO + Krd/n/Kra - K_H/Kra*(Y_H - Y_OH)) \
	#- 0.5*((-Y_O - Y_C + Y_SiO - Krd/n/Kra + K_H/Kra*(Y_H - Y_OH))**2 + 4.*Y_SiO*Y_C)**0.5
	# Removing the H -> OH reaction by eye here, should double check math
	YCO_equil = 0.5*(Y_O + Y_C - Y_SiO + Krd/n/Kra) \
	- 0.5*((-Y_O - Y_C + Y_SiO - Krd/n/Kra)**2 + 4.*Y_SiO*Y_C)**0.5
	return YCO_equil


def Y_SiO_equil(Y_Si, Y_O, Y_CO, n, Kra, Krd):
	'''Solution to quadratic eq for equilibrium number fraction of SiO:
	Y_Si and Y_O are the *total* number fractions.'''
	#YSiO_equil = 0.5*(Y_Si + Y_O + Krd/Kra/n) - \
	#	0.5*((Y_Si + Y_O + Krd/Kra/n)**2 - 4.*Y_Si*Y_O)**0.5
	YSiO_equil = 0.5*(Krd/(n*Kra) - Y_CO + Y_O + Y_Si) - \
		0.5*(4*Y_Si*(Y_CO - Y_O) +(Krd + n*Kra*(-Y_CO + Y_O + Y_Si))**2/Kra**2/n**2)**0.5
	return YSiO_equil



def dYCO_dt(Y_CO, Y_C, Y_O, n, Kra, Krd):
	'''Define function for the change in Y_CO
	this is a first order ODE
	Y_C and Y_O are the *free* number fractions here. '''
	#dYdt = n*Kra*Y_C*Y_O - Krd*Y_CO - n*K_H*Y_CO*Y_H
	# Removing the H reaction part here
	dYdt = n*Kra*Y_C*Y_O - Krd*Y_CO
	return dYdt

def dYH_dt(Y_CO,Y_H,n,K_H):
	'''This includes that change in Y_H due to
	CO + H -> C + OH'''
	#dYdt = -n*K_H*Y_CO*Y_H
	# Need to remove reation completely but doing this first:
	dYdt = 0
	return dYdt

def dYC_dt(Y_CO,Y_H,n,K_H):
	'''This only includes the formation of C from 
	CO + H -> C + OH'''
	#dYdt = n*K_H*Y_CO*Y_H
	# Need to remove reation completely but doing this first:
	dYdt = 0
	return dYdt

def dYOH_dt(Y_CO,Y_H,n,K_H):
	'''This is the formation of OH from the reaction
	CO + H -> C + OH'''
	#dYdt = n*K_H*Y_CO*Y_H
	# Need to remove reation completely but doing this first:
	dYdt = 0
	return dYdt

def dYSiO_dt(Y_SiO, Y_Si, Y_O, n, Kra, Krd):
	'''Define function for the change in Y_CO
	this is a first order ODE
	Y_C and Y_O are the *free* number fractions here. 
	Rates are different than for CO'''
	dYdt = n*Kra*Y_Si*Y_O - Krd*Y_SiO
	return dYdt

