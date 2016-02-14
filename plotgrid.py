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

T_cs, n, R_cs = np.genfromtxt("runs/%s/thermo.txt"%par.directory, unpack=True,skip_footer=1)
t, X_CO, X_C_free, X_O_free, X_D,X_D_sml,X_D_lrg, int_flag, adap_flag, sat = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)
t, n_CO, n_C_free, n_O_free, n_D, Ncrit = np.genfromtxt("runs/%s/densities.txt"%par.directory, unpack=True, skip_footer=1)
K_ra, K_th, K_nth, J  = np.genfromtxt("runs/%s/rates.txt"%par.directory, unpack=True, skip_footer=1)
dust, alldust, size  = np.loadtxt("runs/%s/dust.txt"%par.directory, unpack=True)

dust = np.array(dust[:len(t)])
alldust = np.array(alldust[:len(t)])
size = np.array(size[:len(t)])


#--------- Bin the sizes: -------------
bincount = 500
logsize = np.log10(size)
logsize[np.where(logsize < 0)] = 0.

sizebins = np.logspace(0.,np.max(logsize),bincount)
# logarithmically spaced - range is 0,max(size) not log)    
sizehist = np.zeros(bincount)

for i in range(1,len(sizebins)-1):
    binthis = np.where((size < sizebins[i+1]) & (size > sizebins[i-1]))
    sizehist[i] = np.sum(dust[binthis])

#---------------------------------------


plt.figure(figsize=(12,8))
plt.subplot(221)

plt.plot(t,X_C_free,label=r'$X_{\rm C}^{\rm free}$',color='#013ADF')
plt.plot(t,X_O_free,label=r'$X_{\rm O}^{\rm free}$',color='#04B45F')
#plt.plot(t[np.where(int_flag > 0)], Y_CO[np.where(int_flag > 0)], linewidth=15, alpha=0.1, color='red')
#plt.plot(t[np.where(adap_flag > 0)], Y_CO[np.where(adap_flag > 0)], linewidth=15, alpha=0.1, color='blue')
plt.plot(t,X_CO,label=r'$X_{\rm CO}$',color='red')
plt.plot(t,X_D,label=r'$X_{\rm C}^{\rm solid}$',color='purple')
plt.xlabel('$t$')
plt.title(r'${\rm mass \, fraction}$')
plt.ylabel('$X$')
plt.xscale('log')
plt.ylim([-0.1,1])
plt.xlim([par.tmin,par.tmax])
plt.legend(loc=2)

plt.subplot(222)
#plt.plot([par.tmin,par.tmax],[1.,1.],'--',color='grey',label=r'$S=1$')
plt.semilogy(t, sat)
plt.xlabel('$t$')
plt.ylabel('$S$')
plt.xscale('log')
plt.legend()
#plt.ylim([1e-10,1e10])

plt.subplot(223)
plt.plot(t, n)
plt.xlabel('$t$')
plt.xscale('log')
plt.yscale('log')


plt.subplot(224)
plt.plot(sizebins,sizehist,linestyle='none',marker='o')
plt.xlabel(r'$\log(N_{\rm atoms})$')
plt.ylabel(r'$\#$')
plt.title(r'$\rm size \, distribution$')
plt.xscale('log')
plt.yscale('log')
#plt.ylim([1e-1,1e8])
#plt.xlim([1e0,1e8])



plt.tight_layout()
plt.show()

#plt.savefig('runs/%s/plot.png'%par.directory)