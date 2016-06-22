# Reading in data, doing whatever, plotting it
# copy and paste this into a separate script with the
# necessary packages and same rc_params
import pylab as plt
import numpy as np
import parameters as par
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
print 'dir = ', par.directory

plt.rcParams['font.size'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['lines.linewidth'] = 2


#---------------------------------------------------------

'''Mass fraction evolution over time'''

#par.directory = fdsafdsa
t, dt, X_CO, X_C, X_O, X_solidC, \
X_Mg,X_Si,X_SiO,X_Mg2SiO4,X_H,\
adap_flag, satC, satSi = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)

t=t/86400. # days!
#plt.figure(figsize=(15,6))

plt.plot(t,X_C,label=r'$X_{\rm C}^{\rm free}$')
plt.plot(t,X_O,label=r'$X_{\rm O}^{\rm free}$')
plt.plot(t,X_CO,label=r'$X_{\rm CO}$')
plt.plot(t,X_Mg,label=r'$X_{\rm Mg}$')
plt.plot(t,X_SiO,label=r'$X_{\rm SiO}$')
plt.plot(t,X_Mg2SiO4,label=r'$X_{\rm Mg2SiO4}$')
plt.plot(t,X_solidC,label=r'$X_{\rm C}^{\rm solid}$')
plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'${\rm mass \, fraction}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-14,1e0])
plt.yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e0])
plt.legend(loc=3)
plt.xlim([5e-1,t[-1]])
plt.tight_layout()
plt.show()
#plt.savefig('/Users/Shark/Desktop/talk/X_%s_dt%.0e_Cx%i_e%.0e.eps'%(par.traj,par.dt_init,par.Camt,par.epsilon),format='eps')


#---------------------------------------------------------

'''Plotting all trajectories'''

file_PR  = "trajectories/PR_fiducial.dat"
file_floorSH = "trajectories/floorSH.dat"
file_floorPR = "trajectories/floorPR.dat"
file_nofloor = "trajectories/nofloor.dat"

#t_PR,dens_PR,temp_PR = np.genfromtxt(file_PR,unpack=True)
t_PR,dens_PR,temp_PR,Lshock_PR = np.genfromtxt(file_PR,unpack=True)
t_floorSH,dens_floorSH,temp_floorSH,Lshock_floorSH = np.genfromtxt(file_floorSH,unpack=True)
t_floorPR,dens_floorPR,temp_floorPR,Lshock_floorPR = np.genfromtxt(file_floorPR,unpack=True)
t_nofloor,dens_nofloor,temp_nofloor,Lshock_nofloor = np.genfromtxt(file_nofloor,unpack=True)

plt.plot(dens_PR,temp_PR,label=r'$\rm P \& R \, fiducial$')
plt.plot(dens_floorSH,temp_floorSH,label=r'$\rm Shock \, Luminosity \, floor$',alpha=0.7)
plt.plot(dens_floorPR,temp_floorPR,label=r'$\rm P \& R \, floor$',alpha=0.7)
plt.plot(dens_nofloor,temp_nofloor,label=r'$\rm no \, floor$',alpha=0.7)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$n \, \rm [cm^{-3}]$')
plt.ylabel(r'$T \, \rm [K]$')
plt.legend(loc=2)
plt.show()


#---------------------------------------------------------

'''Pllotting how the max grain size changes over time, from size_hist.py'''

plt.plot(times[1:],maxsizeCg[1:],marker='o',color='#DF0101',markersize=10,markeredgecolor='black')
plt.yscale('linear')
plt.ylim([1e-4,3e-4])
plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'$a_{\rm max} \, \rm [\mu m]$')
plt.tight_layout()
