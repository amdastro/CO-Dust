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
dir1 = '%s/Cx%i/dt%.1e/e1e-01'%(par.traj,par.Camt,par.dt_init)
dir2 = '%s/Cx%i/dt%.1e/e1e-02'%(par.traj,par.Camt,par.dt_init)
dir3 = '%s/Cx%i/dt%.1e/e1e-03'%(par.traj,par.Camt,par.dt_init)



t_1, dt_1, X_CO_1, X_C_1, X_O_1, X_solidC_1, \
X_Mg_1,X_Si_1,X_SiO_1,X_Mg2SiO4_1,X_H_1,\
adap_flag_1, satC_1, satSi_1 = np.genfromtxt("runs/%s/fractions.txt"%dir1, unpack=True,skip_footer=1)
T_1, n_1 = np.genfromtxt("runs/%s/thermo.txt"%dir1, unpack=True,skip_footer=1)

t_2, dt_2, X_CO_2, X_C_2, X_O_2, X_solidC_2, \
X_Mg_2,X_Si_2,X_SiO_2,X_Mg2SiO4_2,X_H_2,\
adap_flag_2, satC_2, satSi_2 = np.genfromtxt("runs/%s/fractions.txt"%dir2, unpack=True,skip_footer=1)
T_2, n_2 = np.genfromtxt("runs/%s/thermo.txt"%dir2, unpack=True,skip_footer=1)


t_3, dt_3, X_CO_3, X_C_3, X_O_3, X_solidC_3, \
X_Mg_3,X_Si_3,X_SiO_3,X_Mg2SiO4_3,X_H_3,\
adap_flag_3, satC_3, satSi_3 = np.genfromtxt("runs/%s/fractions.txt"%dir3, unpack=True,skip_footer=1)
T_3, n_3 = np.genfromtxt("runs/%s/thermo.txt"%dir3, unpack=True,skip_footer=1)



t_1=t_1/86400. # days!
t_2=t_2/86400. # days!
t_3=t_3/86400. # days!

#ax = plt.figure(figsize=(14,10))
#plt.plot(t_1,X_C_1,label=r'$X_{\rm C}^{\rm free}$',linewidth=2,color='#848484')
#plt.plot(t_1,X_O_1,label=r'$X_{\rm O}^{\rm free}$',linewidth=2,color='#848484')
#plt.plot(t_1,X_CO_1,label=r'$X_{\rm CO}$',linewidth=2,color='#848484')
#plt.plot(t_1,X_Mg_1,label=r'$X_{\rm Mg}$',linewidth=2,color='#848484')
#plt.plot(t_1,X_SiO_1,label=r'$X_{\rm SiO}$',linewidth=2,color='#848484')
plt.plot(t_1,X_Mg2SiO4_1,label=r'$X_{\rm Mg2SiO4}$',color='#08298A',linestyle=':')
plt.plot(t_1,X_solidC_1,label=r'$X_{\rm C}^{\rm solid}$',color='#8A0808',linestyle=':')
plt.plot(t_2,X_Mg2SiO4_2,label=r'$X_{\rm Mg2SiO4}$',color='#08298A',linestyle='--')
plt.plot(t_2,X_solidC_2,label=r'$X_{\rm C}^{\rm solid}$',color='#8A0808',linestyle='--')
plt.plot(t_3,X_Mg2SiO4_3,label=r'$X_{\rm Mg2SiO4}$',color='#08298A')
plt.plot(t_3,X_solidC_3,label=r'$X_{\rm C}^{\rm solid}$',color='#8A0808')
plt.annotate(r'$X_{\rm solid \, C}$',xy=(5e1,8e-14))
plt.annotate(r'$X_{\rm Mg_2SiO_4}$',xy=(3.1e0,8e-14))
plt.xlabel(r'$t \, \rm [days]$')
plt.ylabel(r'${\rm mass \, fraction}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-14,1e0])
plt.yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e0])
plt.xlim([5e-1,t_1[-1]])
plt.annotate(r'$\epsilon = 10^{-3}$',xy=(2e2,5.4e-10))
plt.annotate(r'$\epsilon = 10^{-2}$',xy=(2e2,1.4e-7))
plt.annotate(r'$\epsilon = 10^{-1}$',xy=(2e2,7e-6))
plt.tight_layout()
#ax.set_rasterized(True)
plt.show()
#plt.savefig('/Users/Shark/Dropbox/NovaDust/paper/figures/X_%s_dt%.0e_Cx%i_e.eps'%(par.traj,par.dt_init,par.Camt),format='eps')
