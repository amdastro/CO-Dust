import pylab as plt
import numpy as np
import parameters as par
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
print 'dir = ', par.directory

plt.rcParams['font.size'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['lines.linewidth'] = 2

t,X_CO,X_C_free,X_O_free,X_D,X_D_sml,X_D_lrg,int_flag,adap_flag,sat  = np.genfromtxt("runs/%s/fractions.txt"%par.directory,unpack=True,skip_footer=1)
dust, alldust, size  = np.loadtxt("runs/%s/dust.txt"%par.directory, unpack=True)
dust = np.array(dust[:len(t)])
alldust = np.array(alldust[:len(t)])
size = np.array(size[:len(t)])

#plt.figure(figsize=(15,6))

plt.plot(t,X_C_free,label=r'$X_{\rm C}^{\rm free}$',color='#013ADF')
plt.plot(t,X_O_free,label=r'$X_{\rm O}^{\rm free}$',color='#04B45F')
#plt.plot(t[np.where(int_flag > 0)], Y_CO[np.where(int_flag > 0)], linewidth=15, alpha=0.1, color='red')
#plt.plot(t[np.where(adap_flag > 0)], Y_CO[np.where(adap_flag > 0)], linewidth=15, alpha=0.1, color='blue')
plt.plot(t,X_CO,label=r'$X_{\rm CO}$',color='red')
plt.plot(t,X_D,label=r'$X_{\rm dust}$',color='purple')
plt.xlabel(r'$t \, \rm [s]$')
plt.ylabel(r'${\rm mass \, fraction}$')
plt.xscale('log')
plt.ylim([-0.1,1.1])
plt.xlim([par.tmin,np.max(t)])
plt.legend(loc=5)


plt.tight_layout()
plt.show()
#plt.savefig('../X_eg.eps',format='eps')