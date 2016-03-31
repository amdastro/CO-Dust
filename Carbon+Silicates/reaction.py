import physical as ph
import numpy as np


def COformation(T):
	'''Formation rate from Lazzati 2015'''
	K = 4.467e-17/((T/4467.)**(-2.08) + (T/4467.)**(-0.22))**0.5
	return K


def COthermal(T):
	'''Destructiom rate by thermal processes'''
	K = 4.467e-17/((T/4467.)**(-2.08) + (T/4467.)**(-0.22))**0.5\
		*(ph.h**2 / (2.*np.pi*ph.muCO*ph.kB*T))**(-1.5) * np.exp(-ph.B_CO/(ph.kB*T))
	return K

def SiOformation(T):
	'''Formation rate from Todini & Ferrara 2001
	same as Cherchneff and Dwek 2009'''
	K = 5.52e-18 * T**(0.31)
	return K


def SiOthermal(T):
	'''Destruction rate by thermal electrons from Cherchneff and Dwek 2009
	E = 98600 K'''
	K = 4.4e-10*np.exp(-98600./T)
	return np.minimum(K,1e300)

#Should I include this?
def SiONdest(T):
	'''SiO + N -> SiN + O
	Rate from Cherchneff & Dwek 2009'''
	K = 3.84e-9*np.exp(-35959./T)
	return np.minimum(K,1e300)



