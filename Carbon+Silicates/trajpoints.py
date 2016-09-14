import pylab as plt
import numpy as np
import parameters as par
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
print 'dir = ', par.directory

plt.rcParams['font.size'] = 25
plt.rcParams['legend.fontsize'] = 25
plt.rcParams['lines.linewidth'] = 4

# Read in directories for different epsilons
dir1 = 'PR_sphere/Cx%i/dt%.1e/e%.0e'%(par.Camt,par.dt_init,par.epsilon)
dir2 = 'floorPR/Cx%i/dt%.1e/e%.0e'%(par.Camt,par.dt_init,par.epsilon)
dir3 = 'floorSH/Cx%i/dt%.1e/e%.0e'%(par.Camt,par.dt_init,par.epsilon)


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


days_1=t_1/86400. # days!
days_2=t_2/86400. # days!
days_3=t_3/86400. # days!


ii_c = np.where(X_solidC_1 > 1e-14)[0][0]
ii_s = np.where(X_Mg2SiO4_1 > 1e-14)[0][0]
jj_c = np.where(X_solidC_2 > 1e-14)[0][0]
jj_s = np.where(X_Mg2SiO4_2 > 1e-14)[0][0]
kk_c = np.where(X_solidC_3 > 1e-14)[0][0]
kk_s = np.where(X_Mg2SiO4_3 > 1e-14)[0][0]



#line1 = plt.plot(n_1,T_1,color='grey',zorder=1)
line2 = plt.plot(n_2,T_2,color='black',alpha=0.7,zorder=1)
#line3 = plt.plot(n_3,T_3,color='grey',linestyle='-.',alpha=0.9,zorder=1)
#seq2 = [5, 3, 5, 3]
#line2[0].set_dashes(seq2)
#seq3 = [11, 4, 11, 4]
#line3[0].set_dashes(seq3)
#plt.scatter(n_1[ii_c],T_1[ii_c],marker='d',s=180,color='#8856a7',\
#	label=r'$\rm C \, grains \, form$',\
#	linewidth=1,edgecolor='black',zorder=2)

#color='#8856a7',\
plt.scatter(n_2[jj_c],T_2[jj_c],marker='d',s=780,color='darkcyan',\
	linewidth=1,edgecolor='black',zorder=2)
#plt.scatter(n_3[kk_c],T_3[kk_c],marker='d',s=180,color='#8856a7',\
#	linewidth=1,edgecolor='black',zorder=2)
#plt.scatter(n_1[ii_s],T_1[ii_s],marker='*',s=350,color='#9ebcda',\
#	label=r'$\rm silicates \, form$',\
#	linewidth=1,edgecolor='black',zorder=2)

#color='#9ebcda',\
plt.scatter(n_2[jj_s],T_2[jj_s],marker='*',s=1650,color='red',\
	linewidth=1,edgecolor='black',zorder=2)
#plt.scatter(n_3[kk_s],T_3[kk_s],marker='*',s=350,color='#9ebcda',\
#	linewidth=1,edgecolor='black',zorder=2)
plt.xscale('log')
plt.yscale('log')
#plt.xlabel(r'$n \, \rm [cm^{-3}]$')
#plt.ylabel(r'$T \, \rm [K]$')
plt.title(r'$\rm Thermodynamic \, \, trajectory$')
plt.xlim([1e7,1e15])
plt.ylim([1e1,1e5])
#plt.legend(loc=2,scatterpoints=1,handlelength=0,frameon=False)
#plt.annotate(r'$\rm shock \, luminosity \, floor$', xy=(1e11,3.3e2),rotation=28)
plt.annotate(r'$\rm WD \, heating \, floor$', xy=(3e10,6.e2),rotation=8.5)
#plt.annotate(r'$\rm uniform \, sphere$', xy=(4.3e10,3.e3),rotation=15.1)

plt.show()




