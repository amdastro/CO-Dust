import numpy as np
import parameters as par
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
print 'reading ', par.directory

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

Y_grains_Sig = np.array(Y_grains_Sig[:len(t)])
Y_grains_Sig_tot = np.array(Y_grains_Sig_tot[:len(t)])
sizes_Sig = np.array(sizes_Sig[:len(t)])