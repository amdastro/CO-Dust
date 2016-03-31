import pylab as plt
import numpy as np
import parameters as par
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
print 'dir = ', par.directory

plt.rcParams['font.size'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['axes.color_cycle']='black'

T_cs, n, delta, R_cs, c_s = np.genfromtxt("runs/%s/thermo.txt"%par.directory, unpack=True,skip_footer=1)
t, X_CO, X_C_free, X_O_free, X_D, int_flag, adap_flag, sat = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)
t, n_CO, n_C_free, n_O_free, n_D, Ncrit = np.genfromtxt("runs/%s/densities.txt"%par.directory, unpack=True, skip_footer=1)
K_ra, K_th, K_nth, J  = np.genfromtxt("runs/%s/rates.txt"%par.directory, unpack=True, skip_footer=1)
dust, alldust, size, allcarbon  = np.loadtxt("runs/%s/dust.txt"%par.directory, unpack=True)

dust = np.array(dust[:len(t)])
alldust = np.array(alldust[:len(t)])
size = np.array(size[:len(t)])
allcarbon = np.array(allcarbon[:len(t)])

plt.figure(figsize=(8,8))
plt.subplot(211)
#plt.plot([par.tmin,par.tmax],[1.,1.],'--',color='grey',label=r'$S=1$')
plt.semilogy(t, J)
plt.xlabel('$t$')
plt.ylabel('$J$')
plt.xscale('log')
plt.legend()
plt.xlim([1e6,1e8])
plt.ylim([1e-10,1e2])

jj = np.where(Ncrit > 0)
plt.subplot(212)
plt.plot(t[jj], Ncrit[jj])
plt.xlabel('$t$')
plt.ylabel(r'$N_{\rm crit}$')
plt.xscale('log')
plt.yscale('log')
plt.xlim([1e6,1e8])

plt.show()