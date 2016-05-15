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
X_Mg,X_Si,X_SiO,X_Mg2SiO4,X_H,\
adap_flag, satC, satSi = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)
#plt.figure(figsize=(15,6))

t=t/86400. # days!

plt.plot(t,X_C,label=r'$X_{\rm C}^{\rm free}$',color='#4477AA')
plt.plot(t,X_O,label=r'$X_{\rm O}^{\rm free}$',color='#44AAAA')
plt.plot(t,X_CO,label=r'$X_{\rm CO}$',color='#44AA77')
#plt.plot(t,X_Si,label=r'$X_{\rm Si}$')
plt.plot(t,X_Mg,label=r'$X_{\rm Mg}$',color='#AAAA44')
#plt.plot(t,X_H,label=r'$X_{\rm H}$',alpha=0.8)
plt.plot(t,X_SiO,label=r'$X_{\rm SiO}$',color='#AA7744')
plt.plot(t,X_Mg2SiO4,label=r'$X_{\rm Mg2SiO4}$',color='#AA4455')
plt.plot(t,X_solidC,label=r'$X_{\rm C}^{\rm solid}$',color='#AA4488')
plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'${\rm mass \, fraction}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-10,1e0])
#plt.legend(loc=3)

plt.xlim([5e-1,t[-1]])

#plt.tight_layout()
plt.show()
#plt.savefig('/Users/Shark/Desktop/X_e%.0e_Cx5.eps'%par.epsilon,format='eps')