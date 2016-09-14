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

# Read in directories for different epsilons
dir1 = '%s/Cx1/dt1.0e+03/e%.0e'%(par.traj,par.epsilon)
dir2 = '%s/CxeqO/dt1.0e+03/e%.0e'%(par.traj,par.epsilon)
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


'''
plt.plot(t_1,X_Mg2SiO4_1,label=r'$X_{\rm Mg2SiO4}$',color='#9ebcda',alpha=1)
plt.plot(t_1,X_solidC_1,label=r'$X_{\rm C}^{\rm solid}$',color='#8856a7',alpha=1)
plt.annotate(r'${\rm Mg_2SiO_4}$',xy=(3.5e0,6e-14))
plt.annotate(r'${\rm C}_{\rm solid}$',xy=(8e1,6e-14))
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-14,1e0])
plt.yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e0])
plt.xlabel(r'$t - t_{\rm shock} \, \rm [days]$')
plt.ylabel(r'${\rm mass \, fraction}$')
plt.xlim([1e0,t_1[-1]])
plt.annotate(r'$\rm C:O = 0.11$', xy=(1.5e1,3e-2))
'''


lineSi = plt.plot(t_3,X_Mg2SiO4_3,label=r'$X_{\rm Mg2SiO4}$',linestyle='--',alpha=0.9,color='black')
lineC = plt.plot(t_3,X_solidC_3,label=r'$X_{\rm C}^{\rm solid}$',linestyle='--',alpha=0.9,color='black')
seq = [15, 5, 15, 5]
lineSi[0].set_dashes(seq)
lineC[0].set_dashes(seq)
plt.annotate(r'${\rm Mg_2SiO_4}$',xy=(2e1,8e-13),rotation=45)
plt.annotate(r'${\rm C}_{\rm solid}$',xy=(8e0,1e-13),rotation=45)
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-14,1e0])
plt.yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e0])
plt.xlabel(r'$t - t_{\rm shock} \, \rm [days]$')
plt.ylabel(r'${\rm mass \, fraction} \, X$')
plt.xlim([1e0,t_3[-1]])
#plt.annotate(r'$\rm C:O = 1.12$', xy=(8e0,3e-2))
plt.title(r'$\rm C\!:\!O = 1.12$')


plt.show()



