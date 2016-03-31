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

t, dt, X_CO, X_C, X_O, X_solidC, \
X_Mg,X_Si,X_SiO,X_Mg2SiO4,\
adap_flag, satC, satSi = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)
#plt.figure(figsize=(15,6))

plt.plot(t,X_C,label=r'$X_{\rm C}^{\rm free}$',alpha=0.8)
plt.plot(t,X_O,label=r'$X_{\rm O}^{\rm free}$',alpha=0.8)
plt.plot(t,X_CO,label=r'$X_{\rm CO}$')
plt.plot(t,X_Si,label=r'$X_{\rm Si}$')
plt.plot(t,X_Mg,label=r'$X_{\rm Mg}$')
plt.plot(t,X_SiO,label=r'$X_{\rm SiO}$',alpha=0.7)
#plt.plot(t,X_Mg2SiO4,label=r'$X_{\rm Mg2SiO4}$',color='magenta')
#plt.plot(t,X_solidC,label=r'$X_{\rm C}^{\rm solid}$',color='purple')
plt.xlabel(r'$t \, \rm [s]$')
plt.ylabel(r'${\rm mass \, fraction}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-4,1e0])
plt.xlim([par.tmin,np.max(t)])
plt.legend(loc=3)


plt.tight_layout()
plt.show()
#plt.savefig('../X_eg.eps',format='eps')