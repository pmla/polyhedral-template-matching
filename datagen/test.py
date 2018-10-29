import math
import numpy as np
import itertools
import scipy.spatial
import diamond_templates
import square_facets
import wb


def get_fcc():
	l = [	[  0.            ,  0.766032346285,  0.766032346285],
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
	return np.array(l)

def clockwise(a, b, c):
	return np.dot(np.cross(a, b), c) < 0

def make_clockwise(points, a, b, c):
	if not clockwise(points[a], points[b], points[c]):
		return (b, a, c)
	else:
		return (a, b, c)

def generate_triangulations(points, colours, sqfacets, eqfacets):

	cs = np.array(list(itertools.product([0, 1], repeat=len(sqfacets))))

	points = np.array([e / np.linalg.norm(e) for e in points])

	triangulations = []
	for ts in cs:
		#facets = list(eqfacets)

		facets = []
		for f in eqfacets:
			unique = np.unique(colours[f])
			if len(unique) > 1:
				facets += [f]
			else:
				m = unique[0]
				a, b, c = f
				facets += [(m, a, b)]
				facets += [(m, b, c)]
				facets += [(m, c, a)]

		for c, sqf in zip(ts, sqfacets):
			indices = [0, 1, 2, 3]
			f0 = sqf[np.roll(indices, 0 + c)[:3]]
			f1 = sqf[np.roll(indices, 2 + c)[:3]]
			facets += [f0, f1]
		facets = np.array([make_clockwise(points, *f) for f in facets])
		triangulations += [facets]
	return np.array(triangulations)

def _facets_to_edges(facets):

	edges = []
	for a, b, c in facets:
		edges += [(a, b)]
		edges += [(b, c)]
		edges += [(c, a)]
	return edges

def go():

	points, colours = diamond_templates.get_diamond_cubic_points()
	#points = get_fcc()

	sqfacets = square_facets.get_square_facets(points)
	squares = points[sqfacets]
	eqfacets = diamond_templates.get_equilateral_facets(points)
	triangulations = generate_triangulations(points, colours, sqfacets, eqfacets)

	vertex_colours = np.array([0] * 4 + [1] * (len(points) - 4))

	'''
	for facets in triangulations:
		edges = _facets_to_edges(facets)
		planar_graphs.plot_points(points, edges = edges, colours=vertex_colours)
	'''

	data = []
	for t in triangulations:

		print np.max(np.bincount(t.reshape(-1))),

		facets = t - np.min(t)
		#labelling = canonical_form(facets)
		(best_code, automorphisms) = wb.weinberg(facets, right_only=True,
								vertex_colours=vertex_colours, edge_colours=None)
		#print '\n', labelling

		##print best_code

		data += [best_code]

		#data += [tuple(labelling)]
		#print labelling
		#return
	print
	print len(data)
	print len(set(data))

def _go():

	points, colours = diamond_templates.get_diamond_hexagonal_points()

	lengths = np.linalg.norm(points[:], axis=1)
	print len(lengths)
	scale = np.mean(lengths)
	print scale
	points = points / scale

	for x, y, z in points:
		print '{ %.12f, %.12f, %.12f },' % (x, y, z)
go()
