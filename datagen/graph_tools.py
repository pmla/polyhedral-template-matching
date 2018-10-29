import math
import numpy as np
import scipy.spatial
import itertools
from numpy.core.umath_tests import inner1d
from numpy import sqrt

def get_ico():
	p = (1 + math.sqrt(5)) / 2
	l = [	[ 0, -1, -p],
		[ 0, -1,  p],
		[ 0,  1, -p],
		[ 0,  1,  p],
		[-1, -p,  0],
		[-1,  p,  0],
		[ 1, -p,  0],
		[ 1,  p,  0],
		[-p,  0, -1],
		[-p,  0,  1],
		[ p,  0, -1],
		[ p,  0,  1],
		[0, 0, 0]]
	return np.array(l) / np.linalg.norm(l[0])

def get_fcc():
	l = [
		[ 1,  1,  0],
		[ 0,  1,  1],
		[ 1,  0,  1],
		[-1, -1,  0],
		[ 0, -1, -1],
		[-1,  0, -1],
		[-1,  1,  0],
		[ 0, -1,  1],
		[-1,  0,  1],
		[ 1, -1,  0],
		[ 0,  1, -1],
		[ 1,  0, -1],
		[0, 0, 0]]
	return np.array(l) / np.linalg.norm(l[0])

def get_bcc():

	l = [
		[ 1,  1,  1],
		[-1,  1,  1],
		[ 1,  1, -1],
		[-1, -1,  1],
		[ 1, -1,  1],
		[-1,  1, -1],
		[-1, -1, -1],
		[ 1, -1, -1],
		[ 2,  0,  0],
		[-2,  0,  0],
		[ 0,  2,  0],
		[ 0, -2,  0],
		[ 0,  0,  2],
		[ 0,  0, -2],

		[0, 0, 0]]
	return np.array(l) / np.mean(np.linalg.norm(l[:-1], axis=1))

def get_hcp():

	l = [
		[    3*sqrt(2),  -3*sqrt(6),           0 ],
		[   -6*sqrt(2),           0,           0 ],
		[   -3*sqrt(2),     sqrt(6),  -4*sqrt(3) ],
		[    3*sqrt(2),     sqrt(6),  -4*sqrt(3) ],
		[            0,  -2*sqrt(6),  -4*sqrt(3) ],
		[   -3*sqrt(2),   3*sqrt(6),           0 ],
		[    3*sqrt(2),   3*sqrt(6),           0 ],
		[    6*sqrt(2),           0,           0 ],
		[   -3*sqrt(2),  -3*sqrt(6),           0 ],
		[            0,  -2*sqrt(6),   4*sqrt(3) ],
		[    3*sqrt(2),     sqrt(6),   4*sqrt(3) ],
		[   -3*sqrt(2),     sqrt(6),   4*sqrt(3) ],

		[0, 0, 0]]

	return np.array(l) / np.linalg.norm(l[0])

def get_sc():
	l = [
		[ 0.,  0., -1.],
		[ 0.,  0.,  1.],
		[ 0., -1.,  0.],
		[ 0.,  1.,  0.],
		[-1.,  0.,  0.],
		[ 1.,  0.,  0.],
		[ 0.,  0.,  0.]
	]
	return np.array(l)

def facets_to_edges(facets):

	edges = []
	for (a, b, c) in facets:
		edges += [tuple(sorted([a, b]))]
		edges += [tuple(sorted([b, c]))]
		edges += [tuple(sorted([c, a]))]
	return np.array(sorted(list(set(edges))))

def add_centre_subtract_mean(points):
	points = np.concatenate((points, [[0, 0, 0]]))
	return points - np.mean(points, axis=0)

def calc_plane_norms(triangles):
	plane_normals = np.cross(triangles[:,0] - triangles[:,2], triangles[:,1] - triangles[:,2])
	plane_normals /= np.linalg.norm(plane_normals, axis=1)[:, np.newaxis]
	return plane_normals

def get_paired_squares(facets, points):

	data = []
	for i in range(len(facets)):
		for j in range(i+1, len(facets)):
			if len(set(facets[i]).intersection(set(facets[j]))) == 2:

				indices = np.array([facets[i], facets[j]])
				norms = calc_plane_norms(points[indices])
				if np.abs(np.dot(norms[0], norms[1])) > 0.9:
					data += [(facets[i], facets[j])]
	return data

def make_facets_clockwise(facets, points):

	triangles = points[facets]
	us = triangles[:,1] - triangles[:,0]
	vs = triangles[:,2] - triangles[:,0]

	crosses = np.cross(us, vs)
	dots = inner1d(triangles[:,0], crosses)
	clockwise = dots < 0

	for i, e in enumerate(clockwise):
		if not e:
			(a, b, c) = facets[i]
			facets[i] = (b, a, c)
	return facets

def get_facet_combinations(points):

	facets = np.sort(scipy.spatial.ConvexHull(points).simplices)
	squares = get_paired_squares(facets, points)

	_squares = [tuple(e) for f in squares for e in f]
	triangle_facets = [tuple(sorted(e)) for e in facets if tuple(e) not in _squares]

	square_edges = []
	for (a, b) in squares:
		common_edge = set(a).intersection(set(b))
		common_edge = tuple(sorted(list(common_edge)))
		square_edges += [common_edge]

	sq_facets = []
	for (a, b) in squares:

		perimeter = set(a) | set(b)
		assert(len(perimeter) == 4)

		common_edge = set(a) & set(b)
		assert(len(common_edge) == 2)

		new_edge = perimeter - common_edge
		assert(len(new_edge) == 2)

		common_edge = tuple(sorted(list(common_edge)))
		new_edge = tuple(sorted(list(new_edge)))

		facet_pair0 = [tuple(sorted(list(perimeter - set([common_edge[i]])))) for i in range(2)]
		facet_pair1 = [tuple(sorted(list(perimeter - set([new_edge[i]])))) for i in range(2)]
		sq_facets += [(facet_pair0, facet_pair1)]

	cs = [list(e) for e in itertools.product(*sq_facets)]

	comb = []
	for i, c in enumerate(cs):

		square_facets = [e for f in c for e in f]
		facets = triangle_facets + square_facets
		facets = make_facets_clockwise(np.array(facets), points)
		comb += [(facets, len(triangle_facets))]
	return comb


ideal_fcc = get_fcc()
ideal_hcp = get_hcp()
ideal_bcc = get_bcc()
ideal_ico = get_ico()
ideal_sc = get_sc()

def print_val(f):
	e = " %.15f" % f
	if f >= 0:
		e = " " + e
	return e

if __name__ == "__main__":

	for structure in [ideal_fcc, ideal_hcp, ideal_bcc, ideal_ico, ideal_sc]:

		structure = np.array(list(structure[-1:]) + list(structure[:-1]))
		print 
		for row in structure:
			print "{", ",".join([print_val(e) for e in row]), "},"

