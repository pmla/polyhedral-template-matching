import math
import numpy as np
import itertools
import scipy.spatial


def get_square_facet_components(points):

	facets = np.sort(scipy.spatial.ConvexHull(points).simplices)

	square = []
	for i, f in enumerate(facets):
		p = points[f]
		v = p[1:] - p[0]
		v = [e / np.linalg.norm(e) for e in v]
		dot = np.dot(*v)
		dot = min(1, max(-1, dot))
		angle = np.arccos(dot)
		if abs(math.degrees(angle) - 60) > 1E-3:
			square += [i]

	return facets[square]

def get_long_edges(points, simplices):

	orientated = []
	for s in simplices:

		best = (0, -1)
		for i in range(3):
			r = np.roll(s, i)
			p = points[r]
			d = np.linalg.norm(p[0] - p[1])
			best = max(best, (d, tuple(r)))
		_, r = best
		orientated += [r]
	return np.array(orientated)

def merge_triangles(orientated):

	d = dict()
	for t in orientated:
		key = tuple(sorted(t[:2]))
		if key not in d:
			d[key] = []
		d[key] += [t]

	l = []
	for k, v in d.items():
		l += [np.unique(v)]
	return np.array(l)

def order_squares(points, squares):

	cs = np.array(list(itertools.permutations(range(4))))

	ordered = []
	for s in squares:
		best = (float("inf"), None)
		for c in cs:
			p = points[s[c]]
			q = points[s[np.roll(c, 1)]]
			ds = np.linalg.norm(p - q, axis=1)
			m = max(np.abs(ds - ds[0]))
			best = min(best, (m, list(c), ds))
		m, c, ds = best
		ordered += [s[c]]
	return np.array(ordered)

def get_square_facets(points):

	simplices = get_square_facet_components(points)
	orientated = get_long_edges(points, simplices)

	squares = merge_triangles(orientated)
	ordered = order_squares(points, squares)
	return ordered

