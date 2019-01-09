import math
import numpy as np
import scipy.spatial
import itertools
from numpy.core.umath_tests import inner1d
from numpy import sqrt

def get_ico():
	l = np.array([
			[-0,  0,  1],
			[ 0, -0, -1],
			[-2 * np.sqrt(1./5 * (5./8 - np.sqrt(5) / 8)),  (1+np.sqrt(5)) / (2 * np.sqrt(5)), -np.sqrt(1./5)],
			[ 2 * np.sqrt(1./5 * (5./8 - np.sqrt(5) / 8)), -(1+np.sqrt(5)) / (2 * np.sqrt(5)),  np.sqrt(1./5)],
			[ 0, -np.sqrt(4./5), -np.sqrt(1./5)],
			[ 0,  np.sqrt(4./5),  np.sqrt(1./5)],
			[ 2 * np.sqrt(1./5 * (5./8 + np.sqrt(5)/8)), -(np.sqrt(5) - 1) / (2 * np.sqrt(5)), -np.sqrt(1./5)],
			[-2 * np.sqrt(1./5 * (5./8 + np.sqrt(5)/8)),  (np.sqrt(5) - 1) / (2 * np.sqrt(5)),  np.sqrt(1./5)],
			[-2 * np.sqrt(1./5 * (5./8 + np.sqrt(5)/8)), -(np.sqrt(5) - 1) / (2 * np.sqrt(5)), -np.sqrt(1./5)],
			[ 2 * np.sqrt(1./5 * (5./8 + np.sqrt(5)/8)),  (np.sqrt(5) - 1) / (2 * np.sqrt(5)),  np.sqrt(1./5)],
			[ 2 * np.sqrt(1./5 * (5./8 - np.sqrt(5)/8)),   (1+np.sqrt(5)) / (2 * np.sqrt(5)),  -np.sqrt(1./5)],
			[-2 * np.sqrt(1./5 * (5./8 - np.sqrt(5)/8)),  -(1+np.sqrt(5)) / (2 * np.sqrt(5)),   np.sqrt(1./5)],
			[ 0,  0,  0],
		])
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

def get_graphene():

	sqrt = np.sqrt
	l = np.array([
	[                   0,  -3./11+6*sqrt(3)/11,                   0 ],
	[  -3*sqrt(3)/22+9./11,  -3*sqrt(3)/11+3./22,                   0 ],
	[  -9./11+3*sqrt(3)/22,  -3*sqrt(3)/11+3./22,                   0 ],
	[  -9./11+3*sqrt(3)/22,  -9./22+9*sqrt(3)/11,                   0 ],
	[  -3*sqrt(3)/22+9./11,  -9./22+9*sqrt(3)/11,                   0 ],
	[ -3*sqrt(3)/11+18./11,                   0,                   0 ],
	[  -3*sqrt(3)/22+9./11,  -9*sqrt(3)/11+9./22,                   0 ],
	[  -9./11+3*sqrt(3)/22,  -9*sqrt(3)/11+9./22,                   0 ],
	[ -18./11+3*sqrt(3)/11,                   0,                   0 ],
	[                   0,                   0,                   0 ],
	])
	return l

def get_fluorite_ca():

	sqrt = np.sqrt

	l = np.array([
			[ 1,  1,  1],
			[ 1,  1, -1],
			[ 1, -1,  1],
			[ 1, -1, -1],
			[-1,  1,  1],
			[-1,  1, -1],
			[-1, -1,  1],
			[-1, -1, -1],
			[ 2,  2,  0],
			[ 0,  2,  2],
			[ 2,  0,  2],
			[-2, -2,  0],
			[ 0, -2, -2],
			[-2,  0, -2],
			[-2,  2,  0],
			[ 0, -2,  2],
			[-2,  0,  2],
			[ 2, -2,  0],
			[ 0,  2, -2],
			[ 2,  0, -2],
			[ 0,  0,  0],
		])

	factor = 2 * (sqrt(2) * 3 + sqrt(3)) / 5
	return l / factor


def get_equilateral_facets(points):

	facets = np.sort(scipy.spatial.ConvexHull(points).simplices)
	facets = np.array(sorted([tuple(e) for e in facets]))

	eq = []
	for i, f in enumerate(facets):
		p = points[f]
		v = p[1:] - p[0]
		v = [e / np.linalg.norm(e) for e in v]
		dot = np.dot(*v)
		dot = min(1, max(-1, dot))
		angle = np.arccos(dot)
		if abs(math.degrees(angle) - 60) < 1E-3:
			eq += [i]

	return facets[eq]

def select_inner_facets(points, facets, name):

	if name == 'dhex':
		return facets[[2,0,5,7]]
	elif name == 'dcub':
		return facets[[0,6,5,7]]

	midpoints = np.mean(points[facets], axis=1)
	cs = np.array(list(itertools.combinations(range(len(midpoints)), 4)))

	print midpoints
	print midpoints[[0,6,5,7]]
	asdf

	data = []
	for c in cs:
		p = midpoints[c]
		if all([e[2] < np.max(midpoints[:,2]) * 0.99 for e in p]):
			continue

		dots = [np.dot(p[i], p[j]) for i in range(len(c)) for j in range(i+1, len(c))]
		mdev = max([abs(e - dots[0]) for e in dots])
		data += [(mdev, list(c))]
	for e in sorted(data):
		print e
	data.sort()
	_, c = data[0]
	return facets[c]

def get_inner_points(points, facets, name):

	facets = select_inner_facets(points, facets, name)

	z = [[0,0,0]]
	inner = [np.mean(np.concatenate((z, points[f])), axis=0) for f in facets]
	return np.array(inner)

def colour_vertices(inner, points):

	ds = scipy.spatial.distance.cdist(points, inner)
	colours = np.argmin(ds, axis=1)
	indices = np.argsort(colours)
	points = points[indices]
	colours = np.concatenate(([0, 1, 2, 3], colours[indices]))
	return np.concatenate((inner, points)), colours

def get_diamond_points(points, name):
	facets = get_equilateral_facets(points)
	inner = get_inner_points(points, facets, name)
	return colour_vertices(inner, points)

def get_diamond_cubic_points():
	points, colours = get_diamond_points(ideal_fcc[:-1], 'dcub')
	points /= np.mean(np.linalg.norm(points, axis=1))
	points = np.concatenate((points, [[0, 0, 0]]))
	colours = np.concatenate((colours, [-1]))
	return points, colours

def get_diamond_hexagonal_points():
	points, colours = get_diamond_points(ideal_hcp[:-1], 'dhex')
	points /= np.mean(np.linalg.norm(points, axis=1))
	points = np.concatenate((points, [[0, 0, 0]]))
	colours = np.concatenate((colours, [-1]))
	return points, colours


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

def invert_array(t):

	inverted = np.zeros(len(t)).astype(np.int)
	for i, e in enumerate(t):
		inverted[e] = i
	return np.array(inverted)

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
ideal_dcub, _ = get_diamond_cubic_points()
ideal_dhex, _ = get_diamond_hexagonal_points()
ideal_graphene = get_graphene()
ideal_fluorite_ca = get_fluorite_ca()

import quat_utils
U0 = quat_utils.quaternion_to_rotation_matrix([ sqrt(3)/2, 0, 0, 0.5 ])
U1 = quat_utils.quaternion_to_rotation_matrix([ 0, sqrt(3)/2, 0.5, 0 ])
U2 = quat_utils.quaternion_to_rotation_matrix([ 0, 0.5, sqrt(3)/2, 0 ])

U3 = quat_utils.quaternion_to_rotation_matrix([ sqrt(2)/2, sqrt(2)/2, 0, 0 ])


ideal_hcp_alt = np.dot(ideal_hcp, U0.T)
ideal_graphene_alt = np.dot(ideal_graphene, U0.T)

ideal_dhex_alt1 = np.dot(ideal_dhex, U0.T)
ideal_dhex_alt2 = np.dot(ideal_dhex, U1.T)
ideal_dhex_alt3 = np.dot(ideal_dhex, U2.T)

ideal_dcub_alt1 = np.dot(ideal_dcub, U3.T)


def print_val(f):
	e = " %.15f" % f
	if f >= 0:
		e = " " + e
	return e

if __name__ == "__main__":

	import fundamental_mappings
	import generators

	for structures in [[ideal_fcc], [ideal_hcp, ideal_hcp_alt], [ideal_bcc], [ideal_ico], [ideal_sc], [ideal_dcub, ideal_dcub_alt1], [ideal_dhex, ideal_dhex_alt1, ideal_dhex_alt2, ideal_dhex_alt3], [ideal_graphene, ideal_graphene_alt]][1:2]:

		structures = [np.array(list(structure[-1:]) + list(structure[:-1])) for structure in structures]

		'''
		structure = np.array([
		[0,0,0],
		[    3*sqrt(2),  -3*sqrt(6),           0 ],
		[   -6*sqrt(2),           0,           0 ],
		[   -3*sqrt(2),   3*sqrt(6),           0 ],
		[    3*sqrt(2),   3*sqrt(6),           0 ],
		[    6*sqrt(2),           0,           0 ],
		[   -3*sqrt(2),  -3*sqrt(6),           0 ],
		])
		'''

		#print
		#for row in structure:
		#	print "{", ",".join([print_val(e) for e in row]), "},"

		fundamental_mappings.find(structures, generators.generator_hcp_conventional)
		asdf

		M = np.dot(structure.T, structure)
		mpi_scale = np.trace(M) / 3
		print mpi_scale

		inv = structure / mpi_scale
		pinv = np.linalg.pinv(structure).T
		assert np.linalg.norm(inv - pinv) < 1E-9
		continue

		lines = []
		for p in structure / mpi_scale:
			line = "{ " + ", ".join([("%.15f" % e).rjust(18) for e in p]) + " }"
			#line = line.replace(".000000000000", ".            ")
			#line = line.replace(".500000000000", ".5           ")
			lines += [line]
		print ",\n".join(lines)
		print

