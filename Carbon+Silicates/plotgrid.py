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

T_cs, n = np.genfromtxt("runs/%s/thermo.txt"%par.directory, unpack=True,skip_footer=1)

t, dt, X_CO, X_C, X_O, X_solidC, \
X_Mg,X_Si,X_SiO,X_Mg2SiO4,X_H,\
adap_flag, satC, satSi \
= np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)

t, n_CO, n_C_free, n_O_free, n_solidC,\
n_Mg,n_Si,n_SiO,n_Mg2SiO4,n_H,ncrit_Cg,ncrit_Sig\
 = np.genfromtxt("runs/%s/densities.txt"%par.directory, unpack=True, skip_footer=1)

K_CO_ra, K_CO_th, K_CO_nth, \
K_SiO_ra, K_SiO_th, K_SiO_nth,K_CO_H,J_Cg,J_Sig  \
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
logsizeCg = np.log10(sizes_Cg)
logsizeCg[np.where(logsizeCg < 0)] = 0.

sizebinsCg = np.logspace(0.,np.max(logsizeCg),bincount)
# logarithmically spaced - range is 0,max(size) not log)    
sizehistCg = np.zeros(bincount)

for i in range(1,len(sizebinsCg)-1):
    binthis = np.where((sizes_Cg < sizebinsCg[i+1]) & (sizes_Cg > sizebinsCg[i-1]))
    sizehistCg[i] = np.sum(Y_grains_Cg[binthis])


logsizeSig = np.log10(sizes_Sig)
logsizeSig[np.where(logsizeSig < 0)] = 0.

sizebinsSig = np.logspace(0.,np.max(logsizeSig),bincount)
# logarithmically spaced - range is 0,max(size) not log)    
sizehistSig = np.zeros(bincount)

for i in range(1,len(sizebinsSig)-1):
    binthis = np.where((sizes_Sig < sizebinsSig[i+1]) & (sizes_Sig > sizebinsSig[i-1]))
    sizehistSig[i] = np.sum(Y_grains_Sig[binthis])
'''
#---------------------------------------


plt.figure(figsize=(16,5))
plt.subplot(131)

plt.plot(t,X_C,label=r'$X_{\rm C}^{\rm free}$',alpha=0.8)
plt.plot(t,X_O,label=r'$X_{\rm O}^{\rm free}$',alpha=0.8)
plt.plot(t,X_CO,label=r'$X_{\rm CO}$',alpha=0.8)
#plt.plot(t,X_Si,label=r'$X_{\rm Si}$')
plt.plot(t,X_Mg,label=r'$X_{\rm Mg}$',alpha=0.8)
plt.plot(t,X_H,label=r'$X_{\rm H}$',alpha=0.8)
plt.plot(t,X_SiO,label=r'$X_{\rm SiO}$',alpha=0.7)
#plt.plot(t,X_Mg2SiO4,label=r'$X_{\rm Mg2SiO4}$',linewidth=4,alpha=0.7)
#plt.plot(t,X_solidC,label=r'$X_{\rm C}^{\rm solid}$',linewidth=4,alpha=0.7)
plt.xlabel(r'$t \, \rm [s]$')
plt.ylabel(r'${\rm mass \, fraction}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-9,1e0])
plt.xlim([par.tmin,np.max(t)])
plt.legend(loc=3,fontsize=10)


plt.subplot(132)
#plt.plot([par.tmin,par.tmax],[1.,1.],'--',color='grey',label=r'$S=1$')
#plt.semilogy(t, satC,label=r'$S_C$')
#plt.semilogy(t, satSi,label=r'$S_{Si}$')
plt.plot(t,K_CO_ra*n,label=r'$n K_{\rm ra}$',alpha=0.8)
plt.plot(t,K_CO_th,label=r'$K_{\rm th}$',alpha=0.8)
plt.plot(t,K_CO_nth,label=r'$K_{\rm nth}$',alpha=0.8)
#plt.plot(t,K_CO_H*n,label=r'$n K_{\rm H}$',alpha=0.8)
plt.xlabel('$t$')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=1,fontsize=10)
plt.title(r'$\rm CO \, rates$')
plt.ylim([1e-25,1e15])


plt.subplot(133)
plt.plot(t,K_SiO_ra*n,label=r'$n K_{\rm ra}$',alpha=0.8)
plt.plot(t,K_SiO_th,label=r'$K_{\rm th}$',alpha=0.8)
plt.plot(t,K_SiO_nth,label=r'$K_{\rm nth}$',alpha=0.8)
plt.xlabel('$t$')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=1,fontsize=10)
plt.title(r'$\rm SiO \, rates$')
plt.ylim([1e-25,1e15])


'''
plt.subplot(132)
plt.plot(sizebinsCg,sizehistCg,linestyle='none',marker='o')
plt.xlabel(r'$\log(N_{\rm atoms})$')
plt.ylabel(r'$\#$')
plt.title(r'$\rm size \, distribution \, C$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-20,1e2])

plt.subplot(133)
plt.plot(sizebinsSig,sizehistSig,linestyle='none',marker='o')
plt.xlabel(r'$\log(N_{\rm atoms})$')
plt.ylabel(r'$\#$')
plt.title(r'$\rm size \, distribution \, Mg2SiO4$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-20,1e2])
'''


plt.tight_layout()
plt.show()

#plt.savefig('runs/%s/plot.png'%par.directory)