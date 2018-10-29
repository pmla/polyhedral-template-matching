import numpy as np

def walk(best_code, prev, cur, degree, m, common, vertex_colours, edge_colours):

	m = np.copy(m)
	num_edges = np.sum(m)
	num_vertices = m.shape[0]

	index = [-1 for i in range(num_vertices)]

	index[prev] = 0
	code = [0]
	n = 1

	for it in range(1, num_edges):

		m[prev,cur] = 0

		if index[cur] == -1:
			next = common[prev,cur]
		elif m[cur,prev] == 1:
			next = prev
		else:
			next = common[prev,cur]
			while m[cur,next] == 0:
				next = common[next,cur]

		if index[cur] == -1:
			index[cur] = n
			n += 1

		#if index[cur] < best_code[it]:
		#	return ([-1], None)

		code += [index[cur]]
		if vertex_colours is not None:
			code += [vertex_colours[cur]]
		if edge_colours is not None:
			code += [edge_colours[(prev, cur)]]
		prev = cur
		cur = next

	return (code, index)

def weinberg(facets, right_only=True, vertex_colours=None, edge_colours=None):

	edges = []
	for (a, b, c) in facets:
		edges += [(a, b)]
		edges += [(b, c)]
		edges += [(c, a)]
	degree = np.bincount(np.array(edges).reshape(-1))
	n = len(degree)

	m = np.zeros((n, n)).astype(np.int8)
	for (a, b) in edges:
		m[a,b] = 1
		m[b,a] = 1

	common_r = np.zeros((n, n)).astype(np.int8)
	common_l = np.zeros((n, n)).astype(np.int8)
	for (a, b, c) in facets:
		common_r[a,b] = c
		common_r[b,c] = a
		common_r[c,a] = b

		common_l[b,a] = c
		common_l[c,b] = a
		common_l[a,c] = b

	max_edge_degree = max([(degree[a], degree[b]) for (a, b) in edges])

	automorphisms = []
	best_code = [-1] * len(edges)

	sides = ['right', 'left']
	if right_only:
		sides = ['right']

	for side in sides:
		common = [common_l, common_r][side == 'right']
		for (a, b) in edges:
			if 1 and (degree[a], degree[b]) == max_edge_degree:
				(code, index) = walk(best_code, a, b, degree, m, common, vertex_colours, edge_colours)

				if code > best_code:
					best_code = code
					automorphisms = [index]
				elif code == best_code:
					automorphisms += [index]

	return (tuple(best_code), automorphisms)

