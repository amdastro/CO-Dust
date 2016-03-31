import pylab as plt
import numpy as np
import parameters as par
plt.rcParams['font.size'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['lines.linewidth'] = 2

print 'dir = ', par.directory

T_cs, n, delta, R_cs, c_s = np.genfromtxt("runs/%s/thermo.txt"%par.directory, unpack=True,skip_footer=1)
t, X_CO, X_C_free, X_O_free, X_D, int_flag, adap_flag, sat = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)
t, n_CO, n_C_free, n_O_free, n_D = np.genfromtxt("runs/%s/densities.txt"%par.directory, unpack=True, skip_footer=1)
K_ra, K_th, K_nth  = np.genfromtxt("runs/%s/rates.txt"%par.directory, unpack=True, skip_footer=1)
dust, alldust, size, allcarbon  = np.loadtxt("runs/%s/dust.txt"%par.directory, unpack=True)

dust = np.array(dust[:len(t)])
alldust = np.array(alldust[:len(t)])
size = np.array(size[:len(t)])
allcarbon = np.array(allcarbon[:len(t)])



plt.figure(figsize=(12,8))

plt.plot(t,n_C_free,label=r'$n_{\rm C}^{\rm free}$',color='#013ADF')
plt.plot(t,n_O_free,label=r'$n_{\rm O}^{\rm free}$',color='#04B45F')
#plt.plot(t[np.where(int_flag > 0)], Y_CO[np.where(int_flag > 0)], linewidth=15, alpha=0.1, color='red')
#plt.plot(t[np.where(adap_flag > 0)], Y_CO[np.where(adap_flag > 0)], linewidth=15, alpha=0.1, color='blue')
plt.plot(t,n_CO,label=r'$n_{\rm CO}$',color='red')
plt.plot(t,n_D,label=r'$n_{\rm dust}$',color='purple')
plt.xlabel('$t$')
plt.ylabel(r'${\rm number \, density}$')
plt.xscale('log')
plt.yscale('log')
#plt.ylim([-0.1,1.1])
plt.xlim([par.tmin,par.tmax])
plt.legend()


plt.show()