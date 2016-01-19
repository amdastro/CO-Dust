import numpy as np
import parameters as par
import pylab as plt
import os
import re
plt.rcParams['font.size'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['lines.linewidth'] = 3


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
dust = dict()
size = dict()
grow = dict()
for i in range(0,len(files)):
	dust[i], size[i], grow[i]  = np.loadtxt("runs/%s/%s"%(par.directory,files[i]), unpack=True)
	dust[i] = np.array(dust[i][:len(t)])
	size[i] = np.array(size[i][:len(t)])
	grow[i] = np.array(grow[i][:len(t)])


bincount = 100

logsize = dict()
sizebins = dict()
sizehist = dict()
growbins = dict()
growhist = dict()
for i in size:
	logsize[i] = np.log10(size[i])
	logsize[i][np.where(logsize[i] < 0)] = 0.
	# logarithmically spaced - range is 0,10^13): 
	sizebins[i] = np.logspace(0.,np.max(logsize[i]),bincount)
	# growth rate is linearly spaced: 
	growbins[i] = np.linspace(0.,np.max(grow[i]),bincount)
	sizehist[i] = np.zeros(bincount)
	growhist[i] = np.zeros(bincount)
	for j in range(1,bincount-1):
    		bin_sizes = np.where((size[i] < sizebins[i][j+1]) & (size[i] > sizebins[i][j-1]))
    		sizehist[i][j] = np.sum(dust[i][bin_sizes])
    		bin_growth = np.where((grow[i] < growbins[i][j+1]) & (grow[i] > growbins[i][j-1]))
    		# sum the number of grains at this growth rate
    		growhist[i][j] = np.sum(dust[i][bin_growth])

#plt.figure(figsize=(8,10))
'''
plt.subplot(211)
plt.plot(sizebins[0],sizehist[0]*sizebins[0],color='blue',alpha=0.5,label=r'$t= 3 \times 10^6$')
plt.plot(sizebins[1],sizehist[1]*sizebins[1],color='green',alpha=0.5,label=r'$t= 6 \times 10^6$')
plt.plot(sizebins[2],sizehist[2]*sizebins[2],color='orange',alpha=0.5,label=r'$t= 9 \times 10^6$')
plt.plot(sizebins[3],sizehist[3]*sizebins[3],color='red',alpha=0.5,label=r'$t= 1.2 \times 10^7$')
plt.plot(sizebins[32],sizehist[32]*sizebins[32],color='black',alpha=0.5,label=r'$t= 1 \times 10^8$')
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-4,1e10])
plt.xlabel(r'$N_{\rm atoms}$')
plt.title(r'$\rm mass \, distribution$')
plt.xlim([1e0,1e12])
plt.legend()
'''

#plt.subplot(212)
plt.plot(growbins[0],growhist[0],color='blue',alpha=0.5,label=r'$t= 3 \times 10^6$')
plt.plot(growbins[1],growhist[1],color='green',alpha=0.5,label=r'$t= 6 \times 10^6$')
plt.plot(growbins[2],growhist[2],color='orange',alpha=0.5,label=r'$t= 9 \times 10^6$')
plt.plot(growbins[3],growhist[3],color='red',alpha=0.5,label=r'$t= 1.2 \times 10^7$')
plt.plot(growbins[32],growhist[32],color='black',alpha=0.5,label=r'$t= 1 \times 10^8$')
plt.xscale('log')
plt.title(r'$\rm growth \, rate \, hist$')
plt.xlabel(r'$N_{\rm atoms} / dt$')
plt.legend()
plt.tight_layout()


plt.show()