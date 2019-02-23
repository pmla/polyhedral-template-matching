import math
import numpy as np
import scipy.spatial
import itertools
from numpy import sqrt
import rmsd
import graph_tools


#TODO: trees do not respect chirality!

def go():

	ps = graph_tools.ideal_fcc
	n = len(ps)
	indices = np.roll(range(n), 1)
	ps = ps[indices]
	print ps

	prev = [[[0]]]
	for i in range(1, 4):

		#print prev[-1]
		row = []
		for p in prev[-1]:

			#print "p", p

			keeps = []
			for j in range(n):
				if j in p: continue
				trial = p + [j]

				if not keeps:
					dmin = float("inf")
				else:
					dmin = min([rmsd.kabsch_rmsd(ps[trial], ps[k]) for k in keeps])
				keep = dmin > 1E-2
				#print trial, dmin, keep
				if keep:
					keeps += [trial]
			row += keeps
		prev += [row]
		print len(row)
	for e in prev:
		print e
	print
go()
