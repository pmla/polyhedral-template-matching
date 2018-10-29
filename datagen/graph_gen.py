import numpy as np
from weinberg import weinberg
from graph_tools import *


def crap_calc_rmsd(P, Q):
	delta = P - Q
	return np.mean(inner1d(delta, delta))

def crap_kabsch(P, Q):
	A = np.dot(P.T, Q)

	V, S, W = np.linalg.svd(A)
	U = np.dot(V, W)
	return U

def invert_array(t):

	inverted = np.zeros(len(t)).astype(np.int)
	for i, e in enumerate(t):
		inverted[e] = i
	return np.array(inverted)

def get_unique_graphs(points, shape=None):

	assert(len(points) in [6, 12, 14])
	if not shape:
		raise Exception("no shape specified")

	data = []
	combs = get_facet_combinations(points)
	colour_codes = []

	for combit, (facets, num_triangles) in enumerate(combs):

		edges = facets_to_edges(facets)
		code, automorphisms = weinberg(facets)
		new = not any([code == e[0] for e in data])

		if shape in ['fcc', 'hcp']:
			deltas = points[edges[:,0]] - points[edges[:,1]]
			distances = np.linalg.norm(deltas, axis=1)
			threshold = np.mean([min(distances), max(distances)])
			colours = distances < threshold

			edge_colours = dict()
			for (a, b), colour in zip(edges, colours):
				edge_colours[(a, b)] = colour
				edge_colours[(b, a)] = colour

			(colour_code, _dummy) = weinberg(facets, edge_colours=edge_colours)
			new = colour_code not in colour_codes

		if shape == 'bcc':
			(colour_code, _dummy) = weinberg(facets, vertex_colours=[0] * 8 + [1] * 6)
			new = colour_code not in colour_codes

		print combit, len(combs)
		m = automorphisms[0]
		automorphisms = [np.array(list(invert_array(e)[m]) + [len(e)]) for e in automorphisms]

		if not new:
			continue

		ideal_points = {'fcc': ideal_fcc, 'bcc': ideal_bcc, 'hcp': ideal_hcp, 'ico': ideal_ico, 'sc': ideal_sc}[shape]

		facets = [tuple(e) for e in np.sort(facets)]
		d = {}
		for i, k in enumerate(facets):
			if shape != 'bcc':
				d[k] = int(i < num_triangles)
			else:
				d[k] = sum([e < 8 for e in k])

		clrs = []
		auts = []
		n = 0
		for aut in automorphisms:

			_facets = [tuple(sorted([aut[e] for e in f])) for f in facets]
			for e in _facets:
				assert(e in facets)

			clr = tuple([d[f] for f in _facets])
			if clr not in clrs:
				clrs += [clr]
				auts += [tuple(aut)]
		print len(clrs), '/', len(automorphisms)

		trial_points = ideal_points.copy()[:-1]
		trial_points += random_perturbation[:len(trial_points)]
		trial_points = add_centre_subtract_mean(trial_points)
		trial_points = trial_points / np.mean(np.linalg.norm(trial_points, axis=1))

		rmsds = []
		for automorphism in automorphisms:
			mapped = ideal_points[automorphism]

			U = crap_kabsch(mapped, trial_points)
			rmsd = crap_calc_rmsd(np.dot(mapped, U), trial_points)
			rmsd = rmsd**0.5
			rmsds += [int(round(1E9*rmsd))]

		unique_rmsds = list(set(rmsds))
		d = dict()
		for i, rmsd in enumerate(rmsds):
			if rmsd not in d:
				d[rmsd] = []
			d[rmsd] += [i]
		unique_automorphisms = [automorphisms[min(e)] for e in d.values()]
		unique_automorphisms = sorted([e.tolist() for e in unique_automorphisms])
		unique_automorphisms = [np.array(e) for e in unique_automorphisms]
		print "num unique:", len(unique_automorphisms)
		assert(len(unique_automorphisms) == len(clrs))

		for e in unique_automorphisms:
			assert tuple(e.tolist()) in auts

		unique_automorphisms = [np.array([0] + list((e[:-1]+1))) for e in unique_automorphisms]	#reorder automorphisms

		data += [(code, edges, facets, unique_automorphisms)]
		if shape in ['fcc', 'hcp', 'bcc']:
			colour_codes += [colour_code]
	return data

_length = 0.02
random_perturbation = np.random.uniform(-1, 1, (15, 3))
for i, e in enumerate(random_perturbation):
	random_perturbation[i] = _length * e / np.linalg.norm(e)

def print_graph_data(_id, num_vertices, edges, facets, automorphisms, automorphism_index, num_automorphisms):

	lid = ".id = %d" % _id
	l4 = "//.automorphisms = {" + ",\n//".join(["{" + ",".join([str(f) for f in list(e)]) + "}" for e in automorphisms]) + "}"
	l5 = ".facets = {%s}" % (','.join(["{%d,%d,%d}" % (a, b, c) for (a, b, c) in facets]))
	l6 = ".automorphism_index = %d" % automorphism_index
	l7 = ".num_automorphisms = %d" % num_automorphisms
	return "{" + ",\n".join([lid, l6, l7, l4, l5]) + "}"

def dat_to_string(dat, name, auts):

	l = ["graph_t graphs_%s[NUM_%s_GRAPHS] = {" % (name, name.upper())]

	for i, (Gref, edges, facets, automorphisms) in enumerate(dat):
		if type(Gref) == tuple:
			num_vertices = max(Gref) + 1
		else:
			num_vertices = len(Gref.degree())

		if len(automorphisms) == 1:
			s = print_graph_data(i, num_vertices, edges, facets, automorphisms, 0, len(automorphisms)) + ','
		else:
			s = print_graph_data(i, num_vertices, edges, facets, automorphisms, len(auts), len(automorphisms)) + ','
			for e in automorphisms:
				auts += [e]
		l += [s]

	l += ['};']
	return "\n\n".join(l), auts

import hashlib
def md5(fname):
	hash_md5 = hashlib.md5()
	with open(fname, "rb") as f:
		for chunk in iter(lambda: f.read(4096), b""):
			hash_md5.update(chunk)
	return hash_md5.hexdigest()

def go():
	structures = [	('sc', ideal_sc),
			('ico', ideal_ico),
			('fcc', ideal_fcc),
			('hcp', ideal_hcp),
			('bcc', ideal_bcc)	][:]

	auts = [range(17)]
	dump = []
	n = 0
	for (name, points) in structures:

		dat = get_unique_graphs(points[:-1], shape=name)
		sizes = sorted([len(automorphisms) for (G, edges, facets, automorphisms) in dat], reverse=True)
		print len(dat), sum(sizes), sizes
		n += sum(sizes)

		s, auts = dat_to_string(dat, name, auts)
		dump += [s]
	print "sum:", n

	for i, e in enumerate(auts):
		e = list(e)
		e += [-1] * (len(auts[0]) - len(e))
		auts[i] = e
	print auts
	auts = np.array(auts)
	print auts.shape

	s = "int8_t automorphisms[%d][17] = {\n\t" % (len(auts)) + "\n\t".join(['{ ' + ", ".join([str(e).rjust(2) for e in aut]) + '},' for aut in auts]) + "\n};\n\n"

	output_string = s + "\n\n".join(dump) + '\n'
	open('dumped_data.txt', 'w').write(output_string)
	checksum = md5('dumped_data.txt')
	print checksum
	assert checksum.lower() == 'c5257bc0f0de736c4632b24be661fe21'
	print "ok"
go()
