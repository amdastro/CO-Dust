import numpy as np
import pylab as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['lines.linewidth'] = 2

#file_PR = "trajectories/4Andrea_PR.dat"
file_PR  = "runs/trajPR.txt"
file_floorSH = "trajectories/4Andrea_floorSH.dat"
file_floorPR = "trajectories/4Andrea_floorPR.dat"
file_nofloor = "trajectories/4Andrea_nofloor.dat"

#t_PR,dens_PR,temp_PR = np.genfromtxt(file_PR,unpack=True)
t_PR,temp_PR,dens_PR = np.genfromtxt(file_PR,unpack=True)
t_floorSH,dens_floorSH,temp_floorSH,Lshock_floorSH = np.genfromtxt(file_floorSH,unpack=True)
t_floorPR,dens_floorPR,temp_floorPR,Lshock_floorPR = np.genfromtxt(file_floorPR,unpack=True)
t_nofloor,dens_nofloor,temp_nofloor,Lshock_nofloor = np.genfromtxt(file_nofloor,unpack=True)

plt.plot(temp_PR,dens_PR,label=r'$\rm PR$')
plt.plot(temp_floorSH,dens_floorSH,label=r'$\rm floorSH$')
plt.plot(temp_floorPR,dens_floorPR,label=r'$\rm floor PR$')
plt.plot(temp_nofloor,dens_nofloor,label=r'$\rm nofloor$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$T$')
plt.ylabel(r'$n$')

plt.show()