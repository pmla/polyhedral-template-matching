import math
import numpy as np
import itertools
import scipy.spatial


def gen_fcc_points():

	points = list(itertools.product([-1, 0, 1], repeat=3))
	points = [e for e in points if sum(np.abs(e)) == 2]
	return np.array(points)

def gen_hcp_points():

	p = [	[    1,     0,     1],
		[    1,     1,     0],
		[    0,     1,     1],
		[   -1,     1,     0],
		[    0,     1,    -1],
		[    1,     0,    -1],
		[    1,    -1,     0],
		[   -1,     0,     1],
		[    0,    -1,     1],
		[-1./3, -4./3, -1./3],
		[-4./3, -1./3, -1./3],
		[-1./3, -1./3, -4./3]]
	return np.array(p)

def get_equilateral_facets(points):

	facets = np.sort(scipy.spatial.ConvexHull(points).simplices)

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

def select_inner_facets(points, facets):

	midpoints = np.mean(points[facets], axis=1)
	cs = np.array(list(itertools.combinations(range(len(midpoints)), 4)))

	data = []
	for c in cs:
		p = midpoints[c]
		dots = [np.dot(p[i], p[j]) for i in range(len(c)) for j in range(i+1, len(c))]
		mdev = max([abs(e - dots[0]) for e in dots])
		data += [(mdev, list(c))]
	data.sort()
	_, c = data[0]
	return facets[c]

def get_inner_points(points, facets):

	facets = select_inner_facets(points, facets)

	z = [[0,0,0]]
	inner = [np.mean(np.concatenate((z, points[f])), axis=0) for f in facets]
	return np.array(inner)

def colour_vertices(inner, points):

	ds = scipy.spatial.distance.cdist(points, inner)
	colours = np.argmin(ds, axis=1)
	indices = np.argsort(colours)
	points = points[indices]
	colours = np.concatenate(([0, 1, 2, 3], colours[indices]))

	#colours = np.array([0] * len(inner) + [1] * len(points))
	return np.concatenate((inner, points)), colours

def get_diamond_points(func):
	points = func()
	facets = get_equilateral_facets(points)
	inner = get_inner_points(points, facets)
	return colour_vertices(inner, points)

def get_diamond_cubic_points():
	return get_diamond_points(gen_fcc_points)

def get_diamond_hexagonal_points():
	return get_diamond_points(gen_hcp_points)

