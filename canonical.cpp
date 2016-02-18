#include <cstdint>
#include <cassert>
#include <cstring>


//#define MAXF 24
#define MAXV 14
#define MAXE 36

static void build_facet_map(int num_facets, int8_t facets[][3], int8_t common[MAXV][MAXV])
{
	memset(common, -1, sizeof(int8_t) * MAXV * MAXV);

	for (int i = 0;i<num_facets;i++)
	{
		int a = facets[i][0];
		int b = facets[i][1];
		int c = facets[i][2];

		assert(common[a][b] == -1);
		assert(common[b][c] == -1);
		assert(common[c][a] == -1);

		common[a][b] = c;
		common[b][c] = a;
		common[c][a] = b;
	}
}

static void weinberg(int num_nodes, int num_edges, int num_facets, int8_t (*facets)[3], int8_t common[MAXV][MAXV], int8_t* best_code, bool& found, int8_t* canonical_labelling, int prev, int cur)
{
	bool m[num_nodes][num_nodes];
	for (int j = 0;j<num_facets;j++)
	{
		int a = facets[j][0];
		int b = facets[j][1];
		int c = facets[j][2];

		m[a][b] = true;
		m[b][a] = true;

		m[b][c] = true;
		m[c][b] = true;

		m[c][a] = true;
		m[a][c] = true;
	}

	int8_t index[num_nodes];
	memset(index, -1, sizeof(int8_t) * num_nodes);
	index[prev] = 0;
	int n = 1;

	bool winning = false;
	int next = -1;
	for (int it=1;it<2*num_edges;it++)
	{
		m[prev][cur] = false;

		if (index[cur] == -1)
		{
			next = common[prev][cur];
		}
		else if (m[cur][prev] == true)
		{
			next = prev;
		}
		else
		{
			next = common[prev][cur];
			while (m[cur][next] == false)
				next = common[next][cur];
		}

		if (index[cur] == -1)
			index[cur] = n++;

		prev = cur;
		cur = next;

		if (!winning && index[cur] < best_code[it])
			return;

		if (winning || index[cur] > best_code[it])
		{
			winning = true;
			best_code[it] = index[cur];
		}
	}

	if (winning)
	{
		memcpy(canonical_labelling, index, sizeof(int8_t) * num_nodes);
		found = true;
	}
}

uint64_t canonical_form(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree, int8_t* canonical_labelling)
{
	bool found = false;
	int8_t common[MAXV][MAXV] = {0};
	int num_edges = 3 * num_facets / 2;
	build_facet_map(num_facets, facets, common);

	int8_t best_degree[3] = {0};
	for (int i = 0;i<num_facets;i++)
	{
		for (int j = 0;j<3;j++)
		{
			int a = facets[i][(j + 0) % 3];
			int b = facets[i][(j + 1) % 3];
			int c = facets[i][(j + 2) % 3];

			int da = degree[a];
			int db = degree[b];
			int dc = degree[c];

			if (da > best_degree[0] || (da == best_degree[0] && db > best_degree[1]) || (da == best_degree[0] && db == best_degree[1] && dc > best_degree[2]))
			{
				best_degree[0] = da;
				best_degree[1] = db;
				best_degree[2] = dc;
			}
		}
	}

	int8_t best_code[2 * MAXE] = {0};
	memset(best_code, -1, sizeof(int8_t) * 2 * num_edges);
	best_code[0] = 0;

	for (int i = 0;i<num_facets;i++)
	{
		for (int j = 0;j<3;j++)
		{
			int a = facets[i][(j + 0) % 3];
			int b = facets[i][(j + 1) % 3];
			int c = facets[i][(j + 2) % 3];

			int da = degree[a];
			int db = degree[b];
			int dc = degree[c];

			if (da == best_degree[0] && db == best_degree[1] && dc == best_degree[2])
				weinberg(num_nodes, num_edges, num_facets, facets, common, best_code, found, canonical_labelling, a, b);
		}
	}

	canonical_labelling[num_nodes] = num_nodes;

	uint64_t hash = 0;
	for (int i = 0;i<2 * num_edges;i++)
	{
		uint64_t e = best_code[i];
		e += i % 8;
		e &= 0xF;
		e <<= (4 * i) % 64;
		hash ^= e;
	}
	return hash;
}

