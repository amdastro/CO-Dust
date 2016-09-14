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

#read this to get the length of the time array
# and the final mass fractions
t,X_solidC,X_Mg2SiO4 \
	= np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True, usecols=(0,5,9), skip_footer=1)

files = []
path = 'runs/%s/'%par.directory
for f in os.listdir(path):
    if os.path.isfile(os.path.join(path,f)) and f.startswith('final_massfile'):
    	files = np.append(files,f)

# This sorts the files based on the integer at the end
# re.split splits the file names with two delimiters - both _ and .
# the [2] element is the number, and it's converted to an integer then sorted
#files = sorted(files, key = lambda x: int(re.split("[_, .]",x)[2]))
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


bincount = 80

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
		#print i,j, binthis, 'Y = ',Y_grains_Sig[i][binthis], 'sum = ', np.sum(Y_grains_Sig[i][binthis])

# So sizehist*size bins is the number (fraction) of atoms that are in a particular grain size
# X = Y * A_i to get mass fraction 
# then multiplying y axis by a total mass will put it in units of mass normalized to some total ejecta mass

a_binsCg  = dict()
a_binsSig = dict()
peaksizeCg = np.zeros(len(files))
peaksizeSig = np.zeros(len(files))
avgsizeCg = np.zeros(len(files))
avgsizeSig = np.zeros(len(files))
for i in range(0,len(files)):
	a_binsCg[i] = (3.*sizebinsCg[i]/(4.*np.pi*1/ph.vC))**(1./3.) *1e4 # microns
	a_binsSig[i] = (3.*sizebinsSig[i]/(4.*np.pi*1/ph.vMg2SiO4))**(1./3.) * 1e4
	if np.max(sizehistCg[i] > 1e-50): 
		peaksizeCg[i] = a_binsCg[i][np.argmax(sizehistCg[i]*sizebinsCg[i])] # peak of mass dist
		avgsizeCg[i] = np.average(a_binsCg[i],weights=sizehistCg[i]*sizebinsCg[i]) # Average by number, not by mass
		#avgsizeCg[i]  = np.sum(a_binsCg[i]*sizehistCg[i]*sizebinsCg[i])/np.sum(sizehistCg[i]) # This is average by mass, but something might be wrong
	if np.max(sizehistSig[i] > 1e-50): 
		peaksizeSig[i] = a_binsSig[i][np.argmax(sizehistSig[i]*sizebinsSig[i])] # peak of mass dist
		avgsizeSig[i] = np.average(a_binsSig[i],weights=sizehistSig[i]*sizebinsSig[i]) # Average by number, not by mass
		#avgsizeSig[i]  = np.sum(a_binsSig[i]*sizehistSig[i]*sizebinsSig[i])/np.sum(sizehistSig[i]) # This is average by mass, but something is wrong


times = np.zeros(len(files))
# Make an array of times of the files
for i in range(0,len(files)):
	firstsplit = files[i].split("_",2)
	#t_s = firstsplit[2].split(".",1)[0]
	#times[i] = int(t_s)/86400. # Days



print ' '
print par.directory
print ' '
print 'avg size Cg = ', avgsizeCg[0]
print 'avg size Sig = ', avgsizeSig[0]
print ' '
print ' '
print 'peak size Cg = ', peaksizeCg[0]
print 'peak size Sig = ', peaksizeSig[0]
print ' '


def surf_area(a):
    surf = 4.*np.pi*a**2
    #return ['%.3f' % z for z in surf]
    return surf



# Mass dist

# the y-axis hist*bins*A_i is in units of mass fraction X 
# let's multiply it by a total mass, say 10^-4, so put it in units of Msun

M_ejecta = 1e-4

a_dMda_Cg  = sizehistCg[0]*sizebinsCg[0]*ph.A_C*M_ejecta
a_dMda_Sig = sizehistSig[0]*sizebinsSig[0]*ph.A_Mg2SiO4*M_ejecta

plt.scatter(a_binsCg[0],sizehistCg[0]*sizebinsCg[0]*ph.A_C*M_ejecta,color='#8856a7',alpha=1,\
	marker='o',s=150,edgecolor='black',label=r'$\rm Carbon$')
plt.scatter(a_binsSig[0],sizehistSig[0]*sizebinsSig[0]*ph.A_Mg2SiO4*M_ejecta,color='#e0ecf4',alpha=1,\
	marker='d',s=160,edgecolor='black',label=r'$\rm Silicates$')
plt.yscale('log')
plt.xscale('log')
plt.yticks([1e-24,1e-20,1e-16,1e-12,1e-8,1e-4,1e0])
plt.ylim([1e-25,1e-3]) # fid and Cx10
#plt.ylim([1e-25,1e0]) # sphere
#plt.xlim([1e-4,1.1e1])  #fid
plt.xlim([1e-4,1.2e2]) #Cx10
#plt.xlim([1e-4,1e-2]) #sphere
#plt.xlabel(r'$a \, \rm [\mu m]$')
#plt.ylabel(r'$a \frac{dM}{da} \, [M_{\odot}]$')
plt.title(r'$\rm MASS$')
plt.legend(loc=2,scatterpoints=1,handlelength=1,frameon=True)
#plt.annotate(r'$\rm C\!:\!O = 0.11$', xy=(1e-2,2e-5))
#plt.annotate(r'$\rm C\!:\!O = 1.12$', xy=(3e-2,2e-5))
#plt.annotate(r'$\rm C\!:\!O = 0.11$', xy=(6e-4,8e-3)) #sphere

plt.show()



# surface area dist
# So now we want to put the y axis in units of cm^2, so see which size grains
# are dominating in surface area
# we want to plot sizehist*sizebins*ph.A_i * 3 * M_ejecta * Msun[g] = cm^2
# Convince yourself that the *3 should be there!

# a in in microns, convert to cm. 1 micron = 1e-4 cm

a_dSda_Cg  = a_dMda_Cg *3*ph.Msun/a_binsCg[0]/1e-4/ph.rhoC
a_dSda_Sig = a_dMda_Sig*3*ph.Msun/a_binsSig[0]/1e-4/ph.rhoMg2SiO4

# But I think you want to plot Y_grains * surface area of that grain size ..... figure this out
# so that would be sizehist * 4 * pi * sizebins**2

'''
plt.scatter(a_binsCg[0],a_dSda_Cg,color='#8856a7',alpha=1,\
	marker='o',s=150,edgecolor='black',label=r'$\rm Carbon$')
plt.scatter(a_binsSig[0],a_dSda_Sig,color='#e0ecf4',alpha=1,\
	marker='d',s=160,edgecolor='black',label=r'$\rm Silicates$')
plt.yscale('log')
plt.xscale('log')
plt.yticks([1e14,1e18,1e22,1e26,1e30,1e34]) # fid
#plt.yticks([1e22,1e26,1e30,1e34,1e38]) # sphere
#plt.ylim([1e20,1e40])  # sphere
plt.ylim([1e14,1e34])  # fid and Cx10
plt.title(r'$\rm SURFACE \, AREA$')
plt.xlim([8e-5,1.1e1]) # fid
#plt.xlim([1e-4,1e2])  # Cx10
#plt.xlim([1e-4,1e-2])  # sphere
plt.xlabel(r'$a \, \rm [\mu m]$')
plt.ylabel(r'$a \frac{d \Sigma}{da} \, \rm [cm^2]$')
plt.legend(loc=2,scatterpoints=1,handlelength=1,frameon=True)
plt.annotate(r'$\rm C\!:\!O = 0.11$', xy=(8e-3,3e32)) # fid
#plt.annotate(r'$\rm C\!:\!O = 1.12$', xy=(2e-2,3e32)) # Cx10
#plt.annotate(r'$\rm C\!:\!O = 0.11$', xy=(6e-4,3e38)) # sphere
plt.show()
'''




#plt.savefig('/Users/Shark/Dropbox/NovaDust/paper/figures/massdist_%s_Cx%i_e1e-04.pdf'%(par.traj,par.Camt))





