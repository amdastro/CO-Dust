import pylab as plt
import numpy as np
import parameters as par
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
print 'dir = ', par.directory

plt.rcParams['font.size'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['lines.linewidth'] = 3

# Read in directories for different epsilons
dir1 = '%s/Cx1/dt%.1e/e%.0e'%(par.traj,par.dt_init,par.epsilon)
dir2 = '%s/CxeqO/dt%.1e/e%.0e'%(par.traj,par.dt_init,par.epsilon)
dir3 = '%s/Cx10/dt%.1e/e%.0e'%(par.traj,par.dt_init,par.epsilon)

t_1, dt_1, X_CO_1, X_C_1, X_O_1, X_solidC_1, \
X_Mg_1,X_Si_1,X_SiO_1,X_Mg2SiO4_1,X_H_1,\
adap_flag_1, satC_1, satSi_1 = np.genfromtxt("runs/%s/fractions.txt"%dir1, unpack=True,skip_footer=1)

t_2, dt_2, X_CO_2, X_C_2, X_O_2, X_solidC_2, \
X_Mg_2,X_Si_2,X_SiO_2,X_Mg2SiO4_2,X_H_2,\
adap_flag_2, satC_2, satSi_2 = np.genfromtxt("runs/%s/fractions.txt"%dir2, unpack=True,skip_footer=1)

t_3, dt_3, X_CO_3, X_C_3, X_O_3, X_solidC_3, \
X_Mg_3,X_Si_3,X_SiO_3,X_Mg2SiO4_3,X_H_3,\
adap_flag_3, satC_3, satSi_3 = np.genfromtxt("runs/%s/fractions.txt"%dir3, unpack=True,skip_footer=1)

t_1=t_1/86400. # days!
t_2=t_2/86400. # days!
t_3=t_3/86400. # days!


plt.figure(figsize=(8,12))

plt.subplot(211)
plt.plot(t_1,X_Mg2SiO4_1,label=r'$X_{\rm Mg2SiO4}$',color='#B40431',alpha=0.9)
plt.plot(t_1,X_solidC_1,label=r'$X_{\rm C}^{\rm solid}$',color='#084B8A',alpha=0.9)
plt.annotate(r'${\rm Mg_2SiO_4}$',xy=(0.1*t_1[np.argmax(X_Mg2SiO4_1)],1.1*np.max(X_Mg2SiO4_1)))
plt.annotate(r'${\rm C}_{\rm solid}$',xy=(0.1*t_1[np.argmax(X_solidC_1)],1.1*np.max(X_solidC_1)))
plt.title(r'$\rm C:O \le 1$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-14,1e0])
plt.yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e0])
plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'${\rm mass \, fraction}$')
plt.xlim([1e0,t_1[-1]])



plt.subplot(212)
plt.plot(t_3,X_Mg2SiO4_3,label=r'$X_{\rm Mg2SiO4}$',color='#B40431',alpha=0.9)
plt.plot(t_3,X_solidC_3,label=r'$X_{\rm C}^{\rm solid}$',color='#084B8A',alpha=0.9)
plt.annotate(r'${\rm Mg_2SiO_4}$',xy=(3e1,5e-5))
plt.annotate(r'${\rm C}_{\rm solid}$',xy=(1e1,2e-2))
plt.title(r'$\rm C:O > 1$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-14,1e0])
plt.yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e0])
plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'${\rm mass \, fraction}$')
plt.xlim([1e0,t_3[-1]])

plt.tight_layout()
plt.show()



