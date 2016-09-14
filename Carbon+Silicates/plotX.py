import pylab as plt
import numpy as np
import parameters as par
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
print 'dir = ', par.directory

plt.rcParams['font.size'] = 25
plt.rcParams['legend.fontsize'] = 25
plt.rcParams['lines.linewidth'] = 3


t, dt, X_CO, X_C, X_O, X_solidC, \
X_Mg,X_Si,X_SiO,X_Mg2SiO4,X_H,\
adap_flag, satC, satSi = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)
T, n = np.genfromtxt("runs/%s/thermo.txt"%par.directory, unpack=True,skip_footer=1)

days=t/86400. # days!


# shaded when n is less than 10^9
ii = np.where(n < 1e9)[0][0] # put back later
dots = 50
xpts = np.zeros(dots)
ypts = np.logspace(-14,0,dots)
xpts[:] = days[ii] # put back later
plt.plot(xpts,ypts,marker='.',markersize=6,linestyle='none',color='#DF0101',alpha=0.9)
#plt.axvline(days[ii],color='#DF0101',alpha=0.9,linestyle=':',linewidth=2)
plt.annotate(r'$n < 10^9 \, \rm cm^{-3}$',xy=(1.08*days[ii],1e-5),rotation=90,color='#DF0101',fontsize=14)
#plt.axvspan(days[810],days[-1],hatch='//',edgecolor='lightgrey',facecolor='none')


plt.plot(days,X_C,label=r'${\rm C}^{\rm free}$',alpha=0.9,color='black')
#plt.annotate(r'${\rm C}$',xy=(6e-1,1e-2)) # sphere
#xy=(2e0,1.1e-2))
plt.annotate(r'${\rm C}$',xy=(1.3e0,2e-7)) # fid
plt.plot(days,X_O,label=r'${\rm O}^{\rm free}$',alpha=0.9,color='black')
#plt.annotate(r'${\rm O}$',xy=(6e-1,1.15e-1)) # sphere 
#xy=(3.4e0,.6e-4))
plt.annotate(r'${\rm O}$',xy=(9e-1,1.1e-1)) #fid
plt.plot(days,X_CO,label=r'${\rm CO}$',alpha=0.9,color='black')
#plt.annotate(r'${\rm CO}$', xy=(5.9e-1,4e-5)) # sphere 
#xy=(3e-1,1.3e-1))
plt.annotate(r'${\rm CO}$', xy=(5.6e-1,3.1e-3)) # fid
plt.plot(days,X_Mg,label=r'${\rm Mg}$',alpha=0.9,color='black')
#plt.annotate(r'${\rm Mg}$', xy=(1.5e0,1e-4)) # sphere
#xy=(4.e-1,1.2e-4))
plt.annotate(r'${\rm Mg}$', xy=(8.3e0,1e-5)) # fid
plt.plot(days,X_SiO,label=r'${\rm SiO}$',alpha=0.9,color='black')
#plt.annotate(r'${\rm SiO}$',xy=(1.5e0,1.4e-3)) # sphere 
#xy=(4e-1,1.4e-3))
plt.annotate(r'${\rm SiO}$', xy=(1.6e0,1.5e-3)) # fid
lineSi = plt.plot(days,X_Mg2SiO4,label=r'${\rm Mg2SiO4}$',linestyle='--',alpha=0.9,color='black')
#plt.annotate(r'${\rm Mg_2SiO_4}$', xy=(1.25e1,1e-12),rotation=60) # sphere
# xy=(2.e1,1e-10))
plt.annotate(r'${\rm Mg_2SiO_4}$', xy=(3.1e0,7e-13),rotation=45) # fid
lineC = plt.plot(days,X_solidC,label=r'${\rm C}^{\rm solid}$',linestyle='--',alpha=0.9,color='black')
#plt.annotate(r'${\rm C_{solid}}$', xy=(5.1e0,2e-13),rotation=63) # sphere
#xy=(2.8e0,1e-10))
plt.annotate(r'${\rm C_{solid}}$',xy=(6.2e1,1e-13),rotation=50) # fid

seq = [15, 5, 15, 5]
lineSi[0].set_dashes(seq)
lineC[0].set_dashes(seq)


plt.xlim([4.8e-1,days[-1]]) 
#plt.xlabel(r'$t - t_{\rm shock} \, \rm [days]$')
#plt.ylabel(r'${\rm mass \, fraction} \, X$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-14,1e0])
plt.yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e0])
#plt.title(r'$\rm C\!:\!O = 0.11$') # fid
#plt.legend(loc=3)

#plt.xlim([2.5e-1,days[-1]])

plt.show()
#plt.savefig('/Users/Shark/Desktop/talk/X_%s_dt%.0e_Cx%i_e%.0e.eps'%(par.traj,par.dt_init,par.Camt,par.epsilon),format='eps')