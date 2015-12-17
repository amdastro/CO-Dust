import numpy as np
import parameters as par

print 'reading ', par.directory

T_cs, n, delta, R_cs, c_s = np.genfromtxt("runs/%s/thermo.txt"%par.directory, unpack=True,skip_footer=1)
t, Y_CO, Y_C_free, Y_O_free, Y_D, int_flag, adap_flag, sat = np.genfromtxt("runs/%s/fractions.txt"%par.directory, unpack=True,skip_footer=1)
t, X_CO, X_C_free, X_O_free, X_D, dNdt = np.genfromtxt("runs/%s/massfrac.txt"%par.directory, unpack=True, skip_footer=1)
K_ra, K_th, K_nth  = np.genfromtxt("runs/%s/rates.txt"%par.directory, unpack=True, skip_footer=1)
dust, alldust, size, allcarbon  = np.loadtxt("runs/%s/dust.txt"%par.directory, unpack=True)

dust = np.array(dust[:len(t)])
alldust = np.array(alldust[:len(t)])
size = np.array(size[:len(t)])
allcarbon = np.array(allcarbon[:len(t)])