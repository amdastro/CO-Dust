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

T_cs, n = np.genfromtxt("runs/%s/thermo.txt"%par.directory, unpack=True,skip_footer=1)

t, dt, X_CO, X_C, X_O, X_solidC, \
X_Mg,X_Si,X_SiO,X_Mg2SiO4,\
adap_flag, satC, satSi \
= np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)

t, n_CO, n_C_free, n_O_free, n_solidC,\
n_Mg,n_Si,n_SiO,n_Mg2SiO4,ncrit_Cg,ncrit_Sig\
 = np.genfromtxt("runs/%s/densities.txt"%par.directory, unpack=True, skip_footer=1)

K_CO_ra, K_CO_th, K_CO_nth, \
K_SiO_ra, K_SiO_th, K_SiO_nth,J_Cg,J_Sig  \
= np.genfromtxt("runs/%s/rates.txt"%par.directory, unpack=True, skip_footer=1)

Y_grains_Cg, Y_grains_Cg_tot, sizes_Cg  \
= np.loadtxt("runs/%s/Cdust.txt"%par.directory, unpack=True)
Y_grains_Sig, Y_grains_Sig_tot, sizes_Sig  \
= np.loadtxt("runs/%s/Sidust.txt"%par.directory, unpack=True)


Y_grains_Cg = np.array(Y_grains_Cg[:len(t)])
Y_grains_Cg_tot = np.array(Y_grains_Cg_tot[:len(t)])
sizes_Cg = np.array(sizes_Cg[:len(t)])

'''
#--------- Bin the sizes: -------------
bincount = 100
logsize = np.log10(size)
logsize[np.where(logsize < 0)] = 0.

sizebins = np.logspace(0.,np.max(logsize),bincount)
# logarithmically spaced - range is 0,max(size) not log)    
sizehist = np.zeros(bincount)

for i in range(1,len(sizebins)-1):
    binthis = np.where((size < sizebins[i+1]) & (size > sizebins[i-1]))
    sizehist[i] = np.sum(dust[binthis])
'''
#---------------------------------------


plt.figure(figsize=(14,4))
plt.subplot(131)

plt.plot(t,X_C,label=r'$X_{\rm C}^{\rm free}$',color='#013ADF')
plt.plot(t,X_O,label=r'$X_{\rm O}^{\rm free}$',color='#04B45F')
plt.plot(t,X_CO,label=r'$X_{\rm CO}$',color='red')
#plt.plot(t,X_Si,label=r'$X_{\rm Si}$',color='orange')
plt.plot(t,X_SiO,label=r'$X_{\rm SiO}$',color='purple')
#plt.plot(t,X_Mg2SiO4,label=r'$X_{\rm Mg2SiO4}$',color='magenta')
#plt.plot(t,X_solidC,label=r'$X_{\rm C}^{\rm solid}$',color='purple')
plt.xlabel('$t$')
plt.title(r'${\rm mass \, fraction}$')
plt.ylabel('$X$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-5,1e0])
plt.xlim([par.tmin,par.tmax])
plt.legend(loc=3,fontsize=10)
plt.subplot(132)
#plt.plot([par.tmin,par.tmax],[1.,1.],'--',color='grey',label=r'$S=1$')
#plt.semilogy(t, satC,label=r'$S_C$')
#plt.semilogy(t, satSi,label=r'$S_{Si}$')
plt.plot(t,K_CO_ra*n,label=r'$n K_{\rm ra}$')
plt.plot(t,K_CO_th,label=r'$K_{\rm th}$')
plt.plot(t,K_CO_nth,label=r'$K_{\rm nth}$')
plt.xlabel('$t$')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=1,fontsize=10)
plt.title(r'$\rm CO \, rates$')
plt.ylim([1e-25,1e15])

plt.subplot(133)
plt.plot(t,K_SiO_ra*n,label=r'$n K_{\rm ra}$')
plt.plot(t,K_SiO_th,label=r'$K_{\rm th}$')
plt.plot(t,K_SiO_nth,label=r'$K_{\rm nth}$')
plt.xlabel('$t$')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=1,fontsize=10)
plt.title(r'$\rm SiO \, rates$')
plt.ylim([1e-25,1e15])

'''
plt.subplot(224)
plt.plot(sizebins,sizehist,linestyle='none',marker='o')
plt.xlabel(r'$\log(N_{\rm atoms})$')
plt.ylabel(r'$\#$')
plt.title(r'$\rm size \, distribution$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-20,1e2])
'''


plt.tight_layout()
plt.show()

#plt.savefig('runs/%s/plot.png'%par.directory)