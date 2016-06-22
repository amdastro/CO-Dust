import pylab as plt
import numpy as np
import parameters as par
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['lines.linewidth'] = 3

# Read in directories for different epsilons
dir1 = 'PR_sphere/Cx%i/dt%.1e/e%.0e'%(par.Camt,par.dt_init,par.epsilon)
dir2 = 'floorPR/Cx%i/dt%.1e/e%.0e'%(par.Camt,par.dt_init,par.epsilon)
dir3 = 'floorSH/Cx%i/dt%.1e/e%.0e'%(par.Camt,par.dt_init,par.epsilon)


t_1 = np.genfromtxt("runs/%s/fractions.txt"%dir1, unpack=True,skip_footer=1,usecols=(0))
T_1, n_1 = np.genfromtxt("runs/%s/thermo.txt"%dir1, unpack=True,skip_footer=1)

t_2 = np.genfromtxt("runs/%s/fractions.txt"%dir2, unpack=True,skip_footer=1,usecols=(0))
T_2, n_2 = np.genfromtxt("runs/%s/thermo.txt"%dir2, unpack=True,skip_footer=1)

t_3 = np.genfromtxt("runs/%s/fractions.txt"%dir3, unpack=True,skip_footer=1,usecols=(0))
T_3, n_3 = np.genfromtxt("runs/%s/thermo.txt"%dir3, unpack=True,skip_footer=1)

t_1=t_1/86400. # days!
t_2=t_2/86400. # days!
t_3=t_3/86400. # days!

'''
plt.plot(t_1,n_1,label=r'$\rm expanding \, sphere$')
plt.plot(t_2,n_2,label=r'$\rm WD \, heating \, floor$',linestyle='--')
plt.plot(t_3,n_3,label=r'$\rm Shock \, luminosity \, floor$',linestyle=':')
plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'$n \, \rm [cm^{-3}]$')
plt.xscale('log')
plt.yscale('log')
plt.xlim([1e-1,t_1[-1]])
plt.ylim([1e5,1e16])
plt.legend(loc=3)
plt.show()
'''
'''
plt.plot(t_1,T_1,label=r'$\rm expanding \, sphere$')
plt.plot(t_2,T_2,label=r'$\rm WD \, heating \, floor$',linestyle='--')
plt.plot(t_3,T_3,label=r'$\rm Shock \, luminosity \, floor$',linestyle=':')
plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'$T \, \rm [K]$')
plt.xscale('log')
plt.yscale('log')
plt.xlim([1e-1,t_1[-1]])
plt.ylim([1e1,1e4])
plt.legend(loc=3)
plt.show()
'''

plt.plot(n_1,T_1,label=r'$\rm expanding \, sphere$')
plt.plot(n_2,T_2,label=r'$\rm WD \, heating \, floor$',linestyle='--')
plt.plot(n_3,T_3,label=r'$\rm Shock \, luminosity \, floor$',linestyle=':')
plt.xlabel(r'$n \, \rm [cm^{-3}]$')
plt.ylabel(r'$T \, \rm [K]$')
plt.xscale('log')
plt.yscale('log')
#plt.xlim([1e-1,t_1[-1]])
#plt.ylim([1e1,1e4])
plt.legend(loc=2)
plt.show()
