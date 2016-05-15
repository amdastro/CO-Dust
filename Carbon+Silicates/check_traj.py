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
n0 = 1.e15
R0 = 1.44e12
v_ej = 3.21e12/ph.daytosec
tmin = 3.*ph.daytosec
tmax = 50.*ph.daytosec
t = tmin

R = R0
T = T0
n = n0

thermofile = open("runs/trajPR.txt", 'w')
thermofile = open("runs/trajPR.txt", 'a')
thermofile.write("# t    T     n\n")

thermofile.write("%.5f %.5e %.5e\n"%(t, T, n))

dt = 1e2

while t < par.tmax:
	t = t+dt
	R, n, T = traj.pontefract(t,dt,R,R0,v_ej,n,n0,T,T0)
	thermofile.write("%.5f %.5e %.5e\n"%(t, T, n))


thermofile.close()



