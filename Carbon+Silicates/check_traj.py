import numpy as np
import pylab as plt
import physical as ph
import reaction as re
import parameters as par
import abundances as ab
import trajectory as traj
import os
import sys


T0 = 6000.
n0 = 1.e13
R0 = 1.44e12
v_ej = 3.21e12/ph.daytosec
tmin = 0.01*ph.daytosec
tmax = 578.70*ph.daytosec

R = R0
T = T0
n = n0

file_floorPR = "trajectories/floorPR.dat"
t_floorPR,dens_floorPR,temp_floorPR,Lshock_floorPR = np.genfromtxt(file_floorPR,unpack=True)

days = np.zeros(15000)
n = np.zeros(15000)
T = np.zeros(15000)

dt = (tmax - tmin)/15000

t = tmin
for i in range(0,15000):
	t = t+dt
	R, n[i], T[i] = traj.pontefract(t,dt,R,R0,v_ej,n,n0,T,T0)
	days[i] = t/ph.daytosec


sec = days*ph.daytosec
arr = np.array([sec,n,T,Lshock_floorPR]).T
np.savetxt("trajectories/PR_sphere.dat",arr,'%.7e',delimiter='   ')


