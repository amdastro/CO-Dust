import numpy as np
import pylab as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['lines.linewidth'] = 3

#file_PR = "trajectories/4Andrea_PR.dat"
file_1  = "trajectories/PR_sphere.dat"
file_2 = "trajectories/floorSH.dat"
file_3 = "trajectories/floorPR.dat"

#t_PR,dens_PR,temp_PR = np.genfromtxt(file_PR,unpack=True)
t_1,n_1,T_1,L_1 = np.genfromtxt(file_1,unpack=True)
t_2,n_2,T_2,L_2 = np.genfromtxt(file_2,unpack=True)
t_3,n_3,T_3,L_3 = np.genfromtxt(file_3,unpack=True)


days_1=t_1/86400. # days!
days_2=t_2/86400. # days!
days_3=t_3/86400. # days!


'''
plt.plot(days_1,n_1,color='grey',label=r'$\rm expanding \, sphere$')
plt.plot(days_2,n_2,color='#0431B4',linestyle='--',label=r'$\rm Shock \, Luminosity \, floor$')
plt.plot(days_3,n_3,color='#B40404',linestyle='-.',label=r'$\rm WD \, heating \, floor$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'$n \, \rm [cm^{-3}]$')
plt.legend(loc=3)
plt.ylim([1e6,1e15])
plt.xlim([1e-1,1e3])
plt.show()
'''
'''
plt.plot(days_1,T_1,color='grey',label=r'$\rm expanding \, sphere$')
plt.plot(days_2,T_2,color='#0431B4',linestyle='--',label=r'$\rm Shock \, Luminosity \, floor$')
plt.plot(days_3,T_3,color='#B40404',linestyle='-.',label=r'$\rm WD \, heating \, floor$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'$T \, \rm [K]$')
plt.legend(loc=3)
plt.ylim([1e1,1e4])
plt.xlim([1e-1,1e3])
plt.show()
'''



plt.plot(n_1,T_1,color='grey',label=r'$\rm expanding \, sphere$')
plt.plot(n_2,T_2,color='#0431B4',linestyle='--',label=r'$\rm Shock \, Luminosity \, floor$')
plt.plot(n_3,T_3,color='#B40404',linestyle='-.',label=r'$\rm WD \, heating \, floor$')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$T \, \rm [K]$')
plt.xlabel(r'$n \, \rm [cm^{-3}]$')
plt.legend(loc=2)
plt.xlim([1e6,1e15])
plt.ylim([1e1,1e5])
plt.show()




