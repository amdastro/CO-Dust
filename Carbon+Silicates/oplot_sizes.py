import numpy as np
import parameters as par
from matplotlib import pyplot as plt
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

# So that the legend doesn't show markers or lines
plt.rcParams['legend.markerscale'] = 0
plt.rcParams['legend.handletextpad'] = 0

# to put shadow around text
import matplotlib.patheffects as PathEffects

# to draw rectangle
import matplotlib.patches as Rectangle


# Read in directories for different epsilons
dir = dict()
dir[0] = '%s/Cx%i/dt%.1e/e1e-02'%(par.traj,par.Camt,par.dt_init)
dir[1] = '%s/Cx%i/dt%.1e/e1e-03'%(par.traj,par.Camt,par.dt_init)
dir[2] = '%s/Cx%i/dt%.1e/e1e-04'%(par.traj,par.Camt,par.dt_init)


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
bincount = 70
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

M_ejecta = 1e-4

# Mass distributions
a_dMda_Cg0  = sizehistCg[0]*sizebinsCg[0]*ph.A_C*M_ejecta
a_dMda_Sig0 = sizehistSig[0]*sizebinsSig[0]*ph.A_Mg2SiO4*M_ejecta
a_dMda_Cg1  = sizehistCg[1]*sizebinsCg[1]*ph.A_C*M_ejecta
a_dMda_Sig1 = sizehistSig[1]*sizebinsSig[1]*ph.A_Mg2SiO4*M_ejecta
a_dMda_Cg2  = sizehistCg[2]*sizebinsCg[2]*ph.A_C*M_ejecta
a_dMda_Sig2 = sizehistSig[2]*sizebinsSig[2]*ph.A_Mg2SiO4*M_ejecta

# Surface density distributions
a_dSda_Cg0  = a_dMda_Cg0 *3*ph.Msun/a_binsCg[0]/1e-4/ph.rhoC
a_dSda_Sig0 = a_dMda_Sig0*3*ph.Msun/a_binsSig[0]/1e-4/ph.rhoMg2SiO4
a_dSda_Cg1  = a_dMda_Cg1 *3*ph.Msun/a_binsCg[1]/1e-4/ph.rhoC
a_dSda_Sig1 = a_dMda_Sig1*3*ph.Msun/a_binsSig[1]/1e-4/ph.rhoMg2SiO4
a_dSda_Cg2  = a_dMda_Cg2 *3*ph.Msun/a_binsCg[2]/1e-4/ph.rhoC
a_dSda_Sig2 = a_dMda_Sig2*3*ph.Msun/a_binsSig[2]/1e-4/ph.rhoMg2SiO4


# MASS DIST COMPARISON PLOT
'''
fig,ax = plt.subplots()
# C, e = 1e-4
ax = plt.scatter(a_binsCg[2],a_dMda_Cg2,color='#8856a7',\
	marker='o',s=150,edgecolor='black',label=r'${\bf e = 10^{-4}}$',alpha=1)
# Si, e = 1e-4
plt.scatter(a_binsSig[2],a_dMda_Sig2,color='#8856a7',\
	marker='d',s=160,edgecolor='black',alpha=1)
# C, e = 1e-3
plt.scatter(a_binsCg[1],a_dMda_Cg1,color='#9ebcda',\
	marker='o',s=150,edgecolor='black',label=r'${\bf e = 10^{-3}}$',alpha=1)
# Si, e = 1e-3
plt.scatter(a_binsSig[1],a_dMda_Sig1,color='#9ebcda',\
	marker='d',s=160,edgecolor='black',alpha=1)
# C, e = 1e-2
plt.scatter(a_binsCg[0],a_dMda_Cg0,color='#e0ecf4',\
	marker='o',s=150,edgecolor='black',label=r'${\bf e = 10^{-2}}$',alpha=1)
# Si, e = 1e-2
plt.scatter(a_binsSig[0],a_dMda_Sig0,color='#e0ecf4',\
	marker='d',s=160,edgecolor='black',alpha=1)


plt.yscale('log')
plt.yscale('log')
plt.xscale('log')
plt.yticks([1e-16,1e-14,1e-12,1e-10,1e-8,1e-6,1e-4])
plt.ylim([1e-17,1e-4])
plt.xlim([1e-4,1.1e1])
plt.xlabel(r'$a \, \rm [\mu m]$')
plt.ylabel(r'$a \frac{dM}{da} \, [M_{\odot}]$')
plt.annotate(r'$\rm Carbon$', xy=(3e-4,1e-11))
plt.annotate(r'$\rm Silicates$', xy=(9e-1,1e-6))
plt.title(r'$\rm MASS$')
#leg = plt.legend(loc=2,numpoints=1,handlelength=0,frameon=False)
#leg.legendPatch.set_path_effects([PathEffects.withSimplePatchShadow()])
txt = plt.annotate(r'${\bf \epsilon = 10^{-4}}$', xy=(2e-4,1e-5),color='#8856a7',\
	path_effects=[PathEffects.withStroke(linewidth=0.8,foreground="black")],\
	fontsize=20)
txt = plt.annotate(r'${\bf \epsilon = 10^{-3}}$', xy=(2e-4,2e-6),color='#9ebcda',\
	path_effects=[PathEffects.withStroke(linewidth=0.8,foreground="black")],\
	fontsize=20)
txt = plt.annotate(r'${\bf \epsilon = 10^{-2}}$', xy=(2e-4,4e-7),color='#e0ecf4',\
	path_effects=[PathEffects.withStroke(linewidth=0.8,foreground="black")],\
	fontsize=20)
#def color_legend_texts(leg):
#    """Color legend texts based on color of corresponding lines"""
#    for line, txt in zip(leg.get_lines(), leg.get_texts()):
#        txt.set_color(line.get_color())  

#color_legend_texts(leg)

# add a rectangle around the text
# bottom left corner = (x,y) and then width, height
rectangle = plt.Rectangle((1.7e-4, 1.9e-7), 0.0015, 0.00005, fc='none')
plt.gca().add_patch(rectangle)

plt.show()
'''






# SURFACE DENS DIST COMPARISON PLOT

fig,ax = plt.subplots()
# Si, e = 1e-4
plt.scatter(a_binsSig[2],a_dSda_Sig2,color='#8856a7',\
	marker='d',s=160,edgecolor='black',alpha=1)
# Si, e = 1e-3
plt.scatter(a_binsSig[1],a_dSda_Sig1,color='#9ebcda',\
	marker='d',s=160,edgecolor='black',alpha=1)
# Si, e = 1e-2
plt.scatter(a_binsSig[0],a_dSda_Sig0,color='#e0ecf4',\
	marker='d',s=160,edgecolor='black',alpha=1)
# C, e = 1e-2
plt.scatter(a_binsCg[0],a_dSda_Cg0,color='#e0ecf4',\
	marker='o',s=150,edgecolor='black',label=r'${\bf e = 10^{-2}}$',alpha=1)
# C, e = 1e-3
plt.scatter(a_binsCg[1],a_dSda_Cg1,color='#9ebcda',\
	marker='o',s=150,edgecolor='black',label=r'${\bf e = 10^{-3}}$',alpha=1)
# C, e = 1e-4
ax = plt.scatter(a_binsCg[2],a_dSda_Cg2,color='#8856a7',\
	marker='o',s=150,edgecolor='black',label=r'${\bf e = 10^{-4}}$',alpha=1)

plt.yscale('log')
plt.yscale('log')
plt.xscale('log')
plt.yticks([1e14,1e18,1e22,1e26,1e30,1e34])
plt.ylim([1e20,1e34])
plt.xlim([1e-4,1.1e1])
plt.xlabel(r'$a \, \rm [\mu m]$')
plt.ylabel(r'$a \frac{d\Sigma}{da} \, [{\rm cm^{2}}]$')
plt.annotate(r'$\rm Carbon$', xy=(3e-4,1e-11))
plt.annotate(r'$\rm Silicates$', xy=(9e-1,1e-6))
plt.title(r'$\rm SURFACE \, AREA$')
#leg = plt.legend(loc=2,numpoints=1,handlelength=0,frameon=False)
#leg.legendPatch.set_path_effects([PathEffects.withSimplePatchShadow()])
txt = plt.annotate(r'${\bf \epsilon = 10^{-4}}$', xy=(2e-4,8e32),color='#8856a7',\
	path_effects=[PathEffects.withStroke(linewidth=0.8,foreground="black")],\
	fontsize=20)
txt = plt.annotate(r'${\bf \epsilon = 10^{-3}}$', xy=(2e-4,1.2e32),color='#9ebcda',\
	path_effects=[PathEffects.withStroke(linewidth=0.8,foreground="black")],\
	fontsize=20)
txt = plt.annotate(r'${\bf \epsilon = 10^{-2}}$', xy=(2e-4,2.2e31),color='#e0ecf4',\
	path_effects=[PathEffects.withStroke(linewidth=0.8,foreground="black")],\
	fontsize=20)
#def color_legend_texts(leg):
#    """Color legend texts based on color of corresponding lines"""
#    for line, txt in zip(leg.get_lines(), leg.get_texts()):
#        txt.set_color(line.get_color())  

#color_legend_texts(leg)

# add a rectangle around the text
# bottom left corner = (x,y) and then width, height
rectangle = plt.Rectangle((1.7e-4, 9e30), .0015, 6e33, fc='none')
plt.gca().add_patch(rectangle)

plt.show()




#plt.savefig('/Users/Shark/Dropbox/NovaDust/paper/figures/Cmass_%s_dt%.0e_Cx%i_e.eps'%(par.traj,par.dt_init,par.Camt),format='eps')

