import numpy as np
import parameters as par
import pylab as plt
import physical as ph
import os
import re
plt.rcParams['font.size'] = 30
plt.rcParams['legend.fontsize'] = 30
plt.rcParams['lines.linewidth'] = 3
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

print 'reading ', par.directory

#read this to get the length of the time arracy
t = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True, usecols=(0),skip_footer=1)

files = []
path = 'runs/%s/'%par.directory
for f in os.listdir(path):
    if os.path.isfile(os.path.join(path,f)) and f.startswith('massfile'):
    	files = np.append(files,f)

# This sorts the files based on the integer at the end
# re.split splits the file names with two delimiters - both _ and .
# the [2] element is the number, and it's converted to an integer then sorted
files = sorted(files, key = lambda x: int(re.split("[_, .]",x)[2]))
# Now we can read in the files in order!

# Y_grains is the number fraction of grains formed in one step, so J[i] * dt[i] / n[i]
Y_grains_Cg  = dict()
Y_grains_Sig = dict()
# sizes stores/updates the size of each group of grains that are formed at each step
sizes_Cg     = dict()
sizes_Sig    = dict()
# dNdt_grow is the growth rate in atoms/second at each step, so need dt to get full growth
dNdt_grow_Cg  = dict()
dNdt_grow_Sig = dict()

# Read in files from each time snapshot and store the dust size things
# The arrays are initially large, so here we cut them off at tmax
for i in range(0,len(files)):
	Y_grains_Cg[i], Y_grains_Sig[i], sizes_Cg[i], sizes_Sig[i], dNdt_grow_Cg[i], dNdt_grow_Sig[i]  = np.loadtxt("runs/%s/%s"%(par.directory,files[i]), unpack=True)
	Y_grains_Cg[i]  = np.array(Y_grains_Cg[i][:len(t)])
	Y_grains_Sig[i] = np.array(Y_grains_Sig[i][:len(t)])
	sizes_Cg[i]  = np.array(sizes_Cg[i][:len(t)])
	sizes_Sig[i] = np.array(sizes_Sig[i][:len(t)])
	dNdt_grow_Cg[i]  = np.array(dNdt_grow_Cg[i][:len(t)])
	dNdt_grow_Sig[i] = np.array(dNdt_grow_Sig[i][:len(t)])


bincount = 100

# log size makes 
logsizeCg   = dict()
logsizeSig  = dict()
sizebinsCg  = dict()
sizebinsSig = dict()
sizehistCg  = dict()
sizehistSig = dict()

for i in range(0,len(files)):

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
		print i,j, binthis, 'Y = ',Y_grains_Sig[i][binthis], 'sum = ', np.sum(Y_grains_Sig[i][binthis])


a_binsCg  = dict()
a_binsSig = dict()
for i in range(0,len(files)):
	a_binsCg[i] = (3.*sizebinsCg[i]/(4.*np.pi*1/ph.vC))**(1./3.) *1e4 # microns
	a_binsSig[i] = (3.*sizebinsSig[i]/(4.*np.pi*1/ph.vMg2SiO4))**(1./3.) * 1e4

#plt.figure(figsize=(12,5))
'''
#plt.subplot(121)
#plt.plot(sizebinsCg[0],sizehistCg[0],alpha=0.5,linestyle='none',marker='o')
plt.plot(a_binsCg[0],sizehistCg[0]*sizebinsCg[0],color='#AA4455',linestyle='none',\
	marker='o',markersize=10,markeredgecolor='black')
#plt.plot(a_binsCg[1],sizehistCg[1]*sizebinsCg[1],alpha=0.7,linestyle='none',marker='o')
plt.yscale('log')
plt.yscale('log')
plt.xscale('log')
plt.yticks([1e-24,1e-20,1e-16,1e-12,1e-8,1e-4,1e0])
plt.ylim([1e-25,1e3])
#plt.ylim([1e-25,1e1])
#plt.xlim([1e0,1e4])
plt.xlabel(r'$a \, \rm [\mu m]$')
plt.ylabel(r'$a \frac{dN}{da}$')
#plt.ylabel(r'$Y_{\rm solid \, C} \, N_{\rm atoms}$')
plt.title(r'$\rm C$')
#plt.title(r'$\rm C \, mass \, distribution$')

'''
#plt.subplot(122)
#plt.plot(sizebinsSig[0],sizehistSig[0],alpha=0.5,linestyle='none',marker='o')
plt.plot(a_binsSig[0],sizehistSig[0]*sizebinsSig[0],color='#AA4488',\
	linestyle='none',marker='o',markersize=10,markeredgecolor='black')
#plt.plot(a_binsSig[1],sizehistSig[1]*sizebinsSig[1],alpha=0.7,linestyle='none',marker='o')
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-25,1e3])
plt.yticks([1e-24,1e-20,1e-16,1e-12,1e-8,1e-4,1e0])
#plt.ylim([1e-25,1e1])
#plt.xlim([1e0,1e4])
plt.xlabel(r'$a \, \rm [\mu m]$')
plt.ylabel(r'$a \frac{dN}{da}$')
#plt.ylabel(r'$Y_{\rm Mg2SiO4} \, N_{\rm molec}$')
plt.title(r'$\rm Mg{}_2SiO{}_4$')
#plt.title(r'$\rm Mg{}_2SiO{}_4 \, mass \, distribution$')



#plt.tight_layout()
plt.show()
plt.savefig('/Users/Shark/Desktop/mass_e%.0e_Cx5.eps'%par.epsilon,format='eps')




