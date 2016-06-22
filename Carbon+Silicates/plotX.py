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


t, dt, X_CO, X_C, X_O, X_solidC, \
X_Mg,X_Si,X_SiO,X_Mg2SiO4,X_H,\
adap_flag, satC, satSi = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)
T, n = np.genfromtxt("runs/%s/thermo.txt"%par.directory, unpack=True,skip_footer=1)

days=t/86400. # days!


# shaded when n is less than 10^9
ii = np.where(n < 1e9)[0][0]
dots = 50
xpts = np.zeros(dots)
ypts = np.logspace(-14,0,dots)
xpts[:] = days[ii]
plt.plot(xpts,ypts,marker='.',markersize=6,linestyle='none',color='#DF0101',alpha=0.7)
#plt.axvline(days[ii],color='#DF0101',alpha=0.7,linestyle=':')
#plt.axvspan(days[810],days[-1],hatch='//',edgecolor='lightgrey',facecolor='none')

plt.plot(days,X_C,label=r'$X_{\rm C}^{\rm free}$',alpha=0.9,color='black')
plt.annotate(r'$X_{\rm C}$',xy=(6e-1,1e-2)) #xy=(1.1e0,2e-7)) #  # xy=(2e0,1.1e-2))#
plt.plot(days,X_O,label=r'$X_{\rm O}^{\rm free}$',alpha=0.9,color='black')
plt.annotate(r'$X_{\rm O}$',xy=(6e-1,1.1e-1)) # xy=(3e-1,1.1e-1)) ## xy=(3.4e0,.6e-4))#
plt.plot(days,X_CO,label=r'$X_{\rm CO}$',alpha=0.9,color='black')
plt.annotate(r'$X_{\rm CO}$',xy=(5.3e-1,5e-5)) # xy=(5.3e-1,2.5e-3)) #xy=(3e-1,1.3e-1))#
plt.plot(days,X_Mg,label=r'$X_{\rm Mg}$',alpha=0.9,color='black')
plt.annotate(r'$X_{\rm Mg}$',xy=(1.5e0,1e-4)) # xy=(8e0,.8e-5)) # xy=(4.e-1,1.2e-4))#
plt.plot(days,X_SiO,label=r'$X_{\rm SiO}$',alpha=0.9,color='black')
plt.annotate(r'$X_{\rm SiO}$',xy=(1.5e0,1.2e-3)) # xy=(1.6e0,1.4e-3)) # xy=(4e-1,1.4e-3))#
plt.plot(days,X_Mg2SiO4,label=r'$X_{\rm Mg2SiO4}$',linestyle='--',alpha=0.9,color='black')
plt.annotate(r'$X_{\rm Mg_2SiO_4}$',xy=(1.3e1,1e-11)) # xy=(3.1e0,1e-11)) # xy=(2.e1,1e-10))# 
plt.plot(days,X_solidC,label=r'$X_{\rm C}^{\rm solid}$',linestyle='--',alpha=0.9,color='black')
plt.annotate(r'$X_{\rm solid \, C}$',xy=(2.8e0,1e-10)) # xy=(1e2,3e-11)) #  xy=(2.8e0,1e-10))#

plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'${\rm mass \, fraction} \, X$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-14,1e0])
plt.yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e0])
#plt.legend(loc=3)

plt.xlim([2.5e-1,days[-1]])

plt.tight_layout()
plt.show()
#plt.savefig('/Users/Shark/Desktop/talk/X_%s_dt%.0e_Cx%i_e%.0e.eps'%(par.traj,par.dt_init,par.Camt,par.epsilon),format='eps')