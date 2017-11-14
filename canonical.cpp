#include <string.h>
#include <cstdint>
#include <cstdbool>
#include <algorithm>
#include "ptm_constants.h"


#define MAXE (3 * PTM_MAX_NBRS - 6)

static bool build_facet_map(int num_facets, int8_t facets[][3], int8_t common[PTM_MAX_NBRS][PTM_MAX_NBRS])
{
	memset(common, -1, sizeof(int8_t) * PTM_MAX_NBRS * PTM_MAX_NBRS);

	for (int i = 0;i<num_facets;i++)
	{
		int a = facets[i][0];
		int b = facets[i][1];
		int c = facets[i][2];

		//assert(common[a][b] == -1);
		//assert(common[b][c] == -1);
		//assert(common[c][a] == -1);
		if (common[a][b] != -1 || common[b][c] != -1 || common[c][a] != -1)
			return false;

		common[a][b] = c;
		common[b][c] = a;
		common[c][a] = b;
	}

	return true;
}

static bool weinberg(int num_nodes, int num_edges, int8_t common[PTM_MAX_NBRS][PTM_MAX_NBRS], int8_t* best_code, int8_t* canonical_labelling, int prev, int cur)
{
	bool m[PTM_MAX_NBRS][PTM_MAX_NBRS];
	memset(m, 0, sizeof(bool) * PTM_MAX_NBRS * PTM_MAX_NBRS);

	int8_t index[PTM_MAX_NBRS];
	memset(index, -1, sizeof(int8_t) * PTM_MAX_NBRS);
	index[prev] = 0;
	int n = 1;

	bool winning = false;
	int next = -1;
	for (int it=1;it<2*num_edges;it++)
	{
		m[prev][cur] = true;

		if (index[cur] == -1)
		{
			next = common[prev][cur];
			index[cur] = n++;
		}
		else if (m[cur][prev] == false)
		{
			next = prev;
		}
		else
		{
			next = common[prev][cur];
			while (m[cur][next] == true)
				next = common[next][cur];
		}

		prev = cur;
		cur = next;

		if (!winning && index[cur] > best_code[it])
			return false;

		if (winning || index[cur] < best_code[it])
		{
			winning = true;
			best_code[it] = index[cur];
		}
	}

	if (winning)
	{
		memcpy(canonical_labelling, index, sizeof(int8_t) * num_nodes);
		return true;
	}

	return false;
}

int canonical_form(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree, int8_t* canonical_labelling, uint64_t* p_hash)
{
	int8_t common[PTM_MAX_NBRS][PTM_MAX_NBRS] = {{0}};
	int num_edges = 3 * num_facets / 2;
	if (!build_facet_map(num_facets, facets, common))
		return -1;

	int8_t best_code[2 * MAXE] = {0};
	memset(best_code, 126, sizeof(int8_t) * 2 * num_edges);
	best_code[0] = 0;

	bool equal = true;
	for (int i = 1;i<num_nodes;i++)
		if (degree[i] != degree[0])
			equal = false;

	if (equal)
	{
		weinberg(num_nodes, num_edges, common, best_code, canonical_labelling, facets[0][0], facets[0][1]);
	}
	else
	{
		uint32_t best_degree = 0;
		for (int i = 0;i<num_facets;i++)
		{
			int a = facets[i][0];
			int b = facets[i][1];
			int c = facets[i][2];

			int da = degree[a];
			int db = degree[b];
			int dc = degree[c];

			best_degree = std::max(best_degree, ((uint32_t)da << 16) | ((uint32_t)db << 8) | ((uint32_t)dc << 0));
			best_degree = std::max(best_degree, ((uint32_t)da << 0) | ((uint32_t)db << 16) | ((uint32_t)dc << 8));
			best_degree = std::max(best_degree, ((uint32_t)da << 8) | ((uint32_t)db << 0) | ((uint32_t)dc << 16));
		}

		for (int i = 0;i<num_facets;i++)
		{
			int a = facets[i][0];
			int b = facets[i][1];
			int c = facets[i][2];

			int da = degree[a];
			int db = degree[b];
			int dc = degree[c];

			if (best_degree == (((uint32_t)da << 16) | ((uint32_t)db << 8) | ((uint32_t)dc << 0)))
				weinberg(num_nodes, num_edges, common, best_code, canonical_labelling, a, b);

			if (best_degree == (((uint32_t)da << 0) | ((uint32_t)db << 16) | ((uint32_t)dc << 8)))
				weinberg(num_nodes, num_edges, common, best_code, canonical_labelling, b, c);

			if (best_degree == (((uint32_t)da << 8) | ((uint32_t)db << 0) | ((uint32_t)dc << 16)))
				weinberg(num_nodes, num_edges, common, best_code, canonical_labelling, c, a);
		}
	}

	for (int i = num_nodes-1;i>=0;i--)
		canonical_labelling[i+1] = canonical_labelling[i] + 1;
	canonical_labelling[0] = 0;

	uint64_t hash = 0;
	for (int i = 0;i<2 * num_edges;i++)
	{
		uint64_t e = best_code[i];
		e += i % 8;
		e &= 0xF;
		e <<= (4 * i) % 64;
		hash ^= e;
	}

	*p_hash = hash;
	return PTM_NO_ERROR;
}

int graph_degree(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree)
{
	memset(degree, 0, sizeof(int8_t) * num_nodes);

	for (int i = 0;i<num_facets;i++)
	{
		int a = facets[i][0];
		int b = facets[i][1];
		int c = facets[i][2];

		degree[a]++;
		degree[b]++;
		degree[c]++;
	}

	int8_t max_degree = 0;
	for (int i = 0;i<num_nodes;i++)
		max_degree = std::max(max_degree, degree[i]);

	return max_degree;
}

