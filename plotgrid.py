import pylab as plt
import numpy as np
import parameters as par

print 'dir = ', par.directory

plt.rcParams['font.size'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['lines.linewidth'] = 2

# why is the last line incomplete???
T_cs, n, delta, R_cs, c_s = np.genfromtxt("runs/%s/thermo.txt"%par.directory, unpack=True,skip_footer=1)
t, Y_CO, Y_C_free, Y_O_free, Y_dust, int_flag, adap_flag, sat = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)
K_ra, K_th, K_nth  = np.genfromtxt("runs/%s/rates.txt"%par.directory, unpack=True, skip_footer=1)
dust, alldust, size, allcarbon  = np.loadtxt("runs/%s/dust.txt"%par.directory, unpack=True)

dust = np.array(dust[:len(t)])
alldust = np.array(alldust[:len(t)])
size = np.array(size[:len(t)])
allcarbon = np.array(allcarbon[:len(t)])

K_rd = K_th + K_nth

plt.figure(figsize=(12,5))
plt.subplot(121)
plt.plot(t,Y_C_free,label=r'$Y_{\rm C}^{\rm free}$',color='#013ADF')
plt.plot(t,Y_O_free,label=r'$Y_{\rm O}^{\rm free}$',color='#04B45F')
plt.plot(t[np.where(int_flag > 0)], Y_CO[np.where(int_flag > 0)], linewidth=15, alpha=0.1, color='red')
plt.plot(t[np.where(adap_flag > 0)], Y_CO[np.where(adap_flag > 0)], linewidth=15, alpha=0.1, color='blue')
plt.plot(t,Y_CO,label=r'$Y_{\rm CO}$',color='red')
plt.plot(t,Y_dust,label=r'$Y_{\rm dust}$',color='purple')
plt.xlabel('$t$')
plt.ylabel(r'${\rm number \, fraction}$')
plt.xscale('log')
plt.ylim([-0.1,1.1])
plt.xlim([par.tmin,par.tmax])
plt.legend(loc=2)

jj = np.where(sat > 1)
plt.subplot(122)
plt.semilogy(t[jj], Y_C_free[jj]*n[jj],label=r'$n_C$')
plt.semilogy(t[jj], sat[jj],label='$S$')
plt.semilogy(t[jj], allcarbon[jj],label=r'$ \rm solid C$')
plt.xlabel('$t$')
#plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-10,1e10])
plt.legend()






# destruction timescale
#t_dest = 1./(K_rd)
# formation timescale
#t_form = 1./(K_ra*n)

#--------------------- PLOT ---------------------#

#plt.plot(t, 1./(K_rd)/t_exp)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('$t$')
#plt.ylabel(r'$t_{\rm eq} /t_{\rm exp}$')
#plt.savefig('t_eq.png')



#plt.subplot(221)
#plt.plot(t, n,color='black')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('$t$')
#plt.ylabel(r'$n_{\rm tot}$')

#plt.subplot(222)
#plt.plot(t, T_cs,color='black')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('$t$')
#plt.ylabel(r'$T_{\rm cs}$')



#plt.plot(t,K_th,label=r'$K_{\rm therm}$',linestyle='--',color='orange')
#plt.plot(t,K_nth,label=r'$K_{\rm nontherm}$',linestyle='--',color='magenta')
#plt.plot(t,(K_rd), label=r'$K_{\rm rd}$',color='red')
#plt.plot(t,(K_ra*n), label=r'$K_{\rm ra} n$',color='blue')
#plt.xlabel('$t$')
#plt.ylabel(r'$K \rm [s^{-1}]$')
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(loc=3)
#plt.ylim([1e-80,1e10])



plt.tight_layout()
plt.show()

#plt.savefig('runs/%s/plot.png'%par.directory)