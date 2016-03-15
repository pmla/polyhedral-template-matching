import ptmmodule
import numpy as np
import struct
import scipy.linalg
import matplotlib.pyplot as plt
import scipy.optimize


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

	while 1:
		for i in range(num_atoms):

			relative_positions = pos[nbrs[i]] - pos[i]
			sqdist = np.linalg.norm(relative_positions, axis=1)
			nearest = np.argsort(sqdist)[:14]
			positions = np.zeros((15,3))
			positions[1:] = relative_positions[nearest]
			(struct, alloy, rmsd, scale, rot, F, F_res, P, U) = ptmmodule.index_structure(positions, calculate_strains=1)
			#(struct, alloy, rmsd, scale, rot) = ptmmodule.index_structure(positions)
			#if struct == 2:
			#	vm, r = calc_fcc_strain(F, F_res, P, U, positions[:13])
			#	strains[i] = vm
			result[i] = struct

	return result
	indices = np.where(result == 2)[0]
	plt.hist(strains[indices], bins=200)
	plt.show()

def go():
	dat_pos = open('../FeCu_positions.dat', 'rb').read()
	dat_nbr = open('../FeCu_nbrs.dat', 'rb').read()

	n = len(dat_pos) / 24
	print "num atoms:", n

	pos = np.array(struct.unpack(n * 3 * "d", dat_pos)).reshape((n, 3))
	nbrs = np.array(struct.unpack(n * 14 * "i", dat_nbr)).reshape((n, 14))

	ptm = run(pos, nbrs)
	print ptm
	print np.bincount(ptm)

if __name__ == "__main__":
	go()
