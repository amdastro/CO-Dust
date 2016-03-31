'''Here are extra snippets of code that you can plug in when needed'''


	
# --------------- Output the mass distribution -----------
# Save the mass distribution of grains at different times
# save size array (number of atoms) and dust array (number of grains)
# convert to mass dist in reading
# Put this within the time loop 

if (sat > 1 and t % 100000 < par.dt_init):
	massarray = np.array([dust_Y, size, dNdt_grow]).T
	print 'saving mass at t = ',t
	np.savetxt("runs/%s/massfile_t_%.0f.txt"%(par.directory,t), massarray, '%.5e', delimiter='   ',\
	header="# dust        size      growth rate")         

# --------------------------------------------------------
