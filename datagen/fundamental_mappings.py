import numpy as np
import itertools
import quat_utils
import scipy.optimize
import scipy.spatial
import rmsd


def invert_array(t):

	inverted = np.zeros(len(t)).astype(np.int)
	for i, e in enumerate(t):
		inverted[e] = i
	return np.array(inverted)

def get_mappings(structure):

	n = len(structure)
	cs = itertools.permutations(range(n), 3)

	keeps = {}

	for c in cs:
		Q = rmsd.kabsch(structure[:3], structure[list(c)]).T
		assert np.linalg.det(Q) > 0.5

		rotated = np.dot(structure, Q.T)
		ds = scipy.spatial.distance.cdist(rotated, structure)

		_, res = scipy.optimize.linear_sum_assignment(ds)
		res = invert_array(res)
		obj = rmsd.rmsd(rotated[res], structure)
		if obj > 1E-3:
			continue

		key = tuple(res)
		if key not in keeps:
			keeps[key] = Q

	return keeps

def find(structure, generators=None):

	print structure
	structure = structure[1:]
	keeps = get_mappings(structure)

	data = []
	for mapping, Q in sorted(keeps.items()):
		q = quat_utils.rotation_matrix_to_quaternion(Q)
		indices = np.where(np.abs(q) < 1E-5)[0]
		q[indices] = 0

		if generators is not None:
			generators = np.array(generators)
			ds0 = np.linalg.norm(q - generators, axis=1)
			ds1 = np.linalg.norm(q + generators, axis=1)
			ds = np.min((ds0, ds1), axis=0)
			assert np.min(ds) < 1E-4
			index = np.argmin(ds)
			data += [(index, mapping, q)]
		else:
			data += [(len(data), mapping, q)]

	indices, _, qs = zip(*sorted(data))
	assert indices == tuple(range(len(indices)))
	for index, mapping, q in sorted(data):

		mapping = [0] + list(np.array(mapping) + 1)
		print "{" + ", ".join([str(e).rjust(2) for e in mapping]) + "},"#, q

	_qs = []
	for q in qs:
		t = [1, 0.1, 0.01, 0.001]
		if np.dot(-q, t) > np.dot(q, t):
			q = -q
		_qs += [(np.dot(q, t), tuple(q))]
	_qs = sorted(_qs, reverse=True)
	qs = np.array([e[1] for e in _qs])
	for q in qs:
		print "{" + ", ".join([("%.14f" % e).rjust(18) for e in q]) + "},"

