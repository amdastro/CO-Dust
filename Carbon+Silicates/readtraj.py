import numpy as np
import pylab as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['lines.linewidth'] = 2

#file_PR = "trajectories/4Andrea_PR.dat"
file_PR  = "trajectories/PR_sphere.dat"
file_floorSH = "trajectories/floorSH.dat"
file_floorPR = "trajectories/floorPR.dat"
file_nofloor = "trajectories/nofloor.dat"

#t_PR,dens_PR,temp_PR = np.genfromtxt(file_PR,unpack=True)
t_PR,dens_PR,temp_PR,Lshock_PR = np.genfromtxt(file_PR,unpack=True)
t_floorSH,dens_floorSH,temp_floorSH,Lshock_floorSH = np.genfromtxt(file_floorSH,unpack=True)
t_floorPR,dens_floorPR,temp_floorPR,Lshock_floorPR = np.genfromtxt(file_floorPR,unpack=True)
t_nofloor,dens_nofloor,temp_nofloor,Lshock_nofloor = np.genfromtxt(file_nofloor,unpack=True)

plt.plot(dens_PR,temp_PR,label=r'$\rm expanding \, sphere$')
plt.plot(dens_floorSH,temp_floorSH,label=r'$\rm Shock \, Luminosity \, floor$',alpha=0.7)
plt.plot(dens_floorPR,temp_floorPR,label=r'$\rm WD \, heating \, floor$',alpha=0.7)
plt.plot(dens_nofloor,temp_nofloor,label=r'$\rm no \, floor$',alpha=0.7)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$n \, \rm [cm^{-3}]$')
plt.ylabel(r'$T \, \rm [K]$')
plt.legend(loc=2)
plt.show()
