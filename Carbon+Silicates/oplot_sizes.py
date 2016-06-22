import numpy as np
import parameters as par
import pylab as plt
import physical as ph
import os
import re
plt.rcParams['font.size'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['lines.linewidth'] = 3
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
print 'reading ', par.directory

# Read in directories for different epsilons
dir = dict()
dir[0] = '%s/Cx%i/dt%.1e/e1e-01'%(par.traj,par.Camt,par.dt_init)
dir[1] = '%s/Cx%i/dt%.1e/e1e-02'%(par.traj,par.Camt,par.dt_init)
dir[2] = '%s/Cx%i/dt%.1e/e1e-03'%(par.traj,par.Camt,par.dt_init)


t_1 = np.genfromtxt("runs/%s/fractions.txt"%dir[0], unpack=True, usecols=(0),skip_footer=1)

# Read in files from each directory
# Y_grains is the number fraction of grains formed in one step, so J[i] * dt[i] / n[i]
Y_grains_Cg  = dict()
Y_grains_Sig = dict()
# sizes stores/updates the size of each group of grains that are formed at each step
sizes_Cg     = dict()
sizes_Sig    = dict()
# dNdt_grow is the growth rate in atoms/second at each step, so need dt to get full growth
dNdt_grow_Cg  = dict()
dNdt_grow_Sig = dict()
for i in range(0,len(dir)):
	Y_grains_Cg[i], Y_grains_Sig[i], sizes_Cg[i], sizes_Sig[i], dNdt_grow_Cg[i], dNdt_grow_Sig[i]  = np.loadtxt("runs/%s/final_massfile_t_final.txt"%(dir[i]), unpack=True)
	Y_grains_Cg[i]  = np.array(Y_grains_Cg[i][:len(t_1)])
	Y_grains_Sig[i] = np.array(Y_grains_Sig[i][:len(t_1)])
	sizes_Cg[i]  = np.array(sizes_Cg[i][:len(t_1)])
	sizes_Sig[i] = np.array(sizes_Sig[i][:len(t_1)])
	dNdt_grow_Cg[i]  = np.array(dNdt_grow_Cg[i][:len(t_1)])
	dNdt_grow_Sig[i] = np.array(dNdt_grow_Sig[i][:len(t_1)])


# Bin the sizes for each grain for each epsilon run
bincount = 100
logsizeCg   = dict()
logsizeSig  = dict()
sizebinsCg  = dict()
sizebinsSig = dict()
sizehistCg  = dict()
sizehistSig = dict()


for i in range(0,len(dir)):

	logsizeCg[i] = np.log10(sizes_Cg[i])
	logsizeCg[i][np.where(logsizeCg[i] < 0)] = 0.

	sizebinsCg[i] = np.logspace(0.,np.max(logsizeCg[i]),bincount)
	# logarithmically spaced - range is 0,max(size)    
	sizehistCg[i] = np.zeros(bincount)

	for j in range(1,bincount-1):
		binthis = np.where((sizes_Cg[i] < sizebinsCg[i][j+1]) & (sizes_Cg[i] > sizebinsCg[i][j-1]))
		sizehistCg[i][j] = np.sum(Y_grains_Cg[i][binthis])
    	# This is sum(J_i*dt_i/n_i) for each size bin, so it's a Y_total of grains for each bin


	logsizeSig[i] = np.log10(sizes_Sig[i])
	logsizeSig[i][np.where(logsizeSig[i] < 0)] = 0.

	sizebinsSig[i] = np.logspace(0.,np.max(logsizeSig[i]),bincount)
	# logarithmically spaced - range is 0,max(size) not log)    
	sizehistSig[i] = np.zeros(bincount)

	for j in range(1,bincount-1):
		binthis = np.where((sizes_Sig[i] < sizebinsSig[i][j+1]) & (sizes_Sig[i] > sizebinsSig[i][j-1]))
		sizehistSig[i][j] = np.sum(Y_grains_Sig[i][binthis])
		#print i,j, binthis, 'Y = ',Y_grains_Sig[i][binthis], 'sum = ', np.sum(Y_grains_Sig[i][binthis])

# Get physical sizes, assuming spherical symmetry
a_binsCg  = dict()
a_binsSig = dict()
for i in range(0,len(dir)):
	a_binsCg[i] = (3.*sizebinsCg[i]/(4.*np.pi*1/ph.vC))**(1./3.) *1e4 # microns
	a_binsSig[i] = (3.*sizebinsSig[i]/(4.*np.pi*1/ph.vMg2SiO4))**(1./3.) * 1e4



plt.plot(a_binsCg[2],sizehistCg[2]*sizebinsCg[2],color='#084B8A',linestyle='none',\
	marker='o',markersize=10,markeredgecolor='black',label=r'$\rm C \, grains$',alpha=0.9)
plt.plot(a_binsSig[2],sizehistSig[2]*sizebinsSig[2],color='#B40431',alpha=0.9,\
	linestyle='none',marker='p',markersize=12,markeredgecolor='black',label=r'$\rm silicates$')
plt.plot(a_binsCg[1],sizehistCg[1]*sizebinsCg[1],color='#0174DF',linestyle='none',alpha=0.7, rasterized=True,\
	marker='o',markersize=10,markeredgecolor='black',label=r'$\epsilon = 10^{-2}$')
plt.plot(a_binsSig[1],sizehistSig[1]*sizebinsSig[1],color='#FA5882',alpha=0.7, rasterized=True,\
	linestyle='none',marker='p',markersize=12,markeredgecolor='black',label=r'$ \epsilon = 10^{-2}$')
plt.plot(a_binsCg[0],sizehistCg[0]*sizebinsCg[0],color='#2ECCFA',linestyle='none',alpha=0.3, rasterized=True,\
	marker='o',markersize=10,markeredgecolor='black',label=r'$ \epsilon = 10^{-1}$')
plt.plot(a_binsSig[0],sizehistSig[0]*sizebinsSig[0],color='#F5A9BC',alpha=0.3, rasterized=True,\
	linestyle='none',marker='p',markersize=12,markeredgecolor='black',label=r'$ \epsilon = 10^{-1}$')
#plt.plot(a_binsCg[1],sizehistCg[1]*sizebinsCg[1],color='#0B610B',linestyle='none',\
#	marker='d',markersize=10,markeredgecolor='black',label=r'$\epsilon = 10^{-2}$')

'''
plt.plot(a_binsSig[0],sizehistSig[0]*sizebinsSig[0],color='#0B2161',linestyle='none',\
	marker='o',markersize=10,markeredgecolor='black',alpha=0.7)
plt.plot(a_binsSig[1],sizehistSig[1]*sizebinsSig[1],color='#013ADF',linestyle='none',\
	marker='o',markersize=10,markeredgecolor='black',alpha=0.5,label=r'$\rm Mg2SiO4$')
plt.plot(a_binsSig[2],sizehistSig[2]*sizebinsSig[2],color='#5882FA',linestyle='none',\
	marker='o',markersize=10,markeredgecolor='black',alpha=0.3)
	'''
plt.yscale('log')
plt.yscale('log')
plt.xscale('log')
plt.yticks([1e-24,1e-20,1e-16,1e-12,1e-8,1e-4,1e0])
plt.ylim([1e-12,1e-3])
plt.xlim([1e-4,1e1])
plt.xlabel(r'$a \, \rm [\mu m]$')
plt.ylabel(r'$a \frac{dN}{da}$')
plt.legend(loc=2,numpoints=1,handlelength=0,ncol=3)
plt.tight_layout()
plt.show()
#plt.savefig('/Users/Shark/Dropbox/NovaDust/paper/figures/Cmass_%s_dt%.0e_Cx%i_e.eps'%(par.traj,par.dt_init,par.Camt),format='eps')

