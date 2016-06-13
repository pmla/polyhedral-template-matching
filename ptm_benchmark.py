import ptmmodule
import numpy as np
import struct
import scipy.linalg
import matplotlib.pyplot as plt
import scipy.optimize

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def plot_points(points, label=False, edges=None):

	fig = plt.figure()#figsize=(16,14))
	fig.set_tight_layout(True)
	ax = fig.add_subplot(111, projection='3d')

	(xs, ys, zs) = zip(*points)
	ax.scatter(xs, ys, zs, c='r')

	if edges:
		for (a, b) in edges:
			p = [points[a], points[b]]
			(xs, ys, zs) = zip(*p)
			plt.plot(xs, ys, zs, c='k')

	if label:
		for p, l in zip(points, label):
			ax.text(p[0], p[1], p[2], l, size=30, color='k')

	(xs, ys, zs) = zip(*points)
	lim = max([abs(e) for e in xs+ys+zs])
	ax.set_xlim(-lim, lim)
	ax.set_ylim(-lim, lim)
	ax.set_zlim(-lim, lim)
	plt.show()



def vonmises(P):

	t = (P[0][0] - P[1][1])**2 + (P[1][1] - P[2][2])**2 + (P[2][2] - P[0][0])**2
	t += 6 * (P[1][2]**2 + P[2][0]**2 + P[0][1]**2)
	return np.sqrt(t / 2)

def calc_fcc_strain(cF, cF_res):

	F = np.array(cF).reshape((3, 3))
	r = np.array(cF_res)
	(U, A) = scipy.linalg.polar(F, side='left')

	dot = np.dot(np.cross(U[0], U[1]), U[2])
	if dot < 0:
		raise Exception("nope")
	return (vonmises(A), sum(r))

def calc_fcc_strain(cF, cF_res, cP, cU, neighbours, mapping):

	points_fcc = [	[  0.            ,  0.766032346285,  0.766032346285],
			[  0.            , -0.766032346285, -0.766032346285],
			[  0.            ,  0.766032346285, -0.766032346285],
			[  0.            , -0.766032346285,  0.766032346285],
			[  0.766032346285,  0.            ,  0.766032346285],
			[ -0.766032346285,  0.            , -0.766032346285],
			[  0.766032346285,  0.            , -0.766032346285],
			[ -0.766032346285,  0.            ,  0.766032346285],
			[  0.766032346285,  0.766032346285,  0.            ],
			[ -0.766032346285, -0.766032346285,  0.            ],
			[  0.766032346285, -0.766032346285,  0.            ],
			[ -0.766032346285,  0.766032346285,  0.            ],
			[  0.            ,  0.            ,  0.            ]	]
	points_fcc = np.array(points_fcc)

	neighbours = np.concatenate((neighbours[1:], neighbours[:1]))
	neighbours -= np.mean(neighbours, axis=0)
	scale = np.mean(np.linalg.norm(neighbours, axis=1))
	for i in range(len(neighbours)):
		neighbours[i] /= scale

	neighbours = neighbours[list(mapping)]
	res = np.linalg.lstsq(points_fcc, neighbours)
	_F = res[0].T

	F = np.array(cF).reshape((3, 3))
	r = np.array(cF_res)
	(U, P) = scipy.linalg.polar(F, side='left')

	dot = np.dot(np.cross(U[0], U[1]), U[2])
	if dot < 0:
		raise Exception("nope")

	assert(np.linalg.norm(F - _F) < 1E-6)

	cU = np.array(cU).reshape((3, 3))
	cP = np.array(cP).reshape((3, 3))
	print cP - P

	print np.linalg.norm(U, axis=1)

	return (vonmises(P), sum(r))

def calc_fcc_strain(cF, cF_res, cP, cU, neighbours):

	F = np.array(cF).reshape((3, 3))
	r = np.array(cF_res)
	(U, P) = scipy.linalg.polar(F, side='left')

	dot = np.dot(np.cross(U[0], U[1]), U[2])
	if dot < 0:
		raise Exception("nope")

	U = np.array(cU).reshape((3, 3))
	P = np.array(cP).reshape((3, 3))

	return (vonmises(P), sum(r))

def run(pos, nbrs):

	num_atoms = len(pos)
	result = np.zeros(num_atoms, int)
	strains = np.zeros(num_atoms).astype(np.double)
	rmsds = np.zeros(num_atoms).astype(np.double)

	if 1:
		for i in range(num_atoms):

			positions = np.concatenate(([pos[i]], pos[nbrs[i][:18]]))
			(struct, alloy, rmsd, scale, rot, F, F_res, P, U, lattice_constant) = ptmmodule.index_structure(positions, calculate_strains=1, topological_ordering=1)
			if 0 and struct == 0 and all([e > 1000 for e in nbrs[i]]):
				print nbrs[i]
				print positions
				positions -= pos[i]
				plot_points(positions[:15])
			#(struct, alloy, rmsd, scale, rot, lattice_constant) = ptmmodule.index_structure(positions)
			#if struct == 2:
			#	vm, r = calc_fcc_strain(F, F_res, P, U, positions[:13])
			#	strains[i] = vm
			result[i] = struct
			rmsds[i] = rmsd

	indices = np.where(rmsds < 0.12)[0]
	kept = np.bincount(result[indices])
	kept[0] += num_atoms - len(indices)
	print "rmsd < 0.12:", kept, sum(kept)

	return result, rmsds
	indices = np.where(result == 2)[0]
	plt.hist(strains[indices], bins=200)
	plt.show()

def go():
	if 1:
		dat_pos = open('test_data/FeCu_positions.dat', 'rb').read()
		dat_nbr = open('test_data/FeCu_nbrs.dat', 'rb').read()
		#dat_pos = open('test_data/fcc_positions.dat', 'rb').read()
		#dat_nbr = open('test_data/fcc_nbrs.dat', 'rb').read()

		n = len(dat_pos) / 24
		print "num atoms:", n

		pos = np.array(struct.unpack(n * 3 * "d", dat_pos)).reshape((n, 3))
		nbrs = np.array(struct.unpack(n * 24 * "i", dat_nbr)).reshape((n, 24))

		n = len(dat_pos) / 24
		print "num atoms:", n

	#pos = np.load('positions0.npy')
	#nbrs = np.load('neighbours0.npy')

	ptm, rmsds = run(pos, nbrs)
	print ptm
	print np.bincount(ptm)

	np.save('structures', ptm)
	np.save('rmsds', rmsds)

if __name__ == "__main__":
	go()
