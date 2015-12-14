import pylab as plt
import numpy as np
import parameters as par

print 'dir = ', par.directory

plt.rcParams['font.size'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['lines.linewidth'] = 2

t,Y_CO,Y_C_free,Y_O_free,Y_D,int_flag,adap_flag,sat  = np.genfromtxt("runs/%s/fractions.txt"%par.directory,unpack=True,skip_footer=2)
dust, alldust, size, allcarbon  = np.loadtxt("runs/%s/dust.txt"%par.directory, unpack=True)
dust = np.array(dust[:len(t)])
alldust = np.array(alldust[:len(t)])
size = np.array(size[:len(t)])
allcarbon = np.array(allcarbon[:len(t)])

plt.figure(figsize=(12,8))

plt.plot(t,Y_C_free,label=r'$Y_{\rm C}^{\rm free}$',color='#013ADF')
plt.plot(t,Y_O_free,label=r'$Y_{\rm O}^{\rm free}$',color='#04B45F')
#plt.plot(t[np.where(int_flag > 0)], Y_CO[np.where(int_flag > 0)], linewidth=15, alpha=0.1, color='red')
#plt.plot(t[np.where(adap_flag > 0)], Y_CO[np.where(adap_flag > 0)], linewidth=15, alpha=0.1, color='blue')
plt.plot(t,Y_CO,label=r'$Y_{\rm CO}$',color='red')
plt.plot(t,Y_D,label=r'$Y_{\rm dust}$',color='purple')
plt.xlabel('$t$')
plt.ylabel(r'${\rm number \, fraction}$')
plt.xscale('log')
plt.ylim([-0.1,1.1])
plt.xlim([par.tmin,par.tmax])
plt.legend(loc=2,fontsize=20)


plt.show()