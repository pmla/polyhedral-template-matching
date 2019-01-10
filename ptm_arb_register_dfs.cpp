#include <cstdint>
#include <cstdbool>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <cassert>
#include <cmath>
#include "ptm_polar.h"
#include "ptm_arb_treefilter.h"


namespace ptm {

static double evaluate(double G1, double G2, double* A)
{
	double E0 = (G1 + G2) / 2;
	double nrmsdsq;
	FastCalcRMSD(A, E0, &nrmsdsq);
	assert(nrmsdsq == nrmsdsq);
	return nrmsdsq;
}

bool check_filter(int level, const treefilter_t* filter, uint8_t* permutation)
{
	if (level >= filter->num_rows)
		return true;

	const uint8_t (*patterns)[3] = &filter->data[filter->row_starts[level]];

	int num_patterns = filter->row_lengths[level];
	for (int i=0;i<num_patterns;i++)
	{
		bool match = true;
		for (int j=0;j<level+1;j++)
			if (patterns[i][j] != permutation[j])
				match = false;
		if (match)
			return true;
	}

	return false;
}

static int recurse(int num_points, int level, double (*P)[3], double (*Q)[3], double max_obj, const treefilter_t* filter,
			uint8_t* permutation, double* p_G1, double* p_G2, double* A,
			uint8_t* best_permutation, double* p_best_obj)
{
	if (level == num_points)
	{
		double obj = evaluate(*p_G1, *p_G2, A);
		if (obj < *p_best_obj)
		{
			*p_best_obj = obj;
			memcpy(best_permutation, permutation, num_points * sizeof(uint8_t));
		}

		return 1;
	}

	if (level >= 2)
	{
		double obj = evaluate(*p_G1, *p_G2, A);
		if (obj > *p_best_obj || obj > max_obj)
			return 1;
	}

	int num_nodes = 0;
	int num_swaps = num_points - level;
	for (int i=0;i<num_swaps;i++)
	{
		int temp = permutation[level];
		permutation[level] = permutation[level + i];
		permutation[level + i] = temp;

		if (check_filter(level, filter, permutation))
		{
			increment_innerproduct(A, level, P, Q, permutation, p_G1, p_G2);

			int num_rnodes = recurse(num_points, level + 1, P, Q, max_obj, filter,
							permutation, p_G1, p_G2, A,
							best_permutation, p_best_obj);
			num_nodes += num_rnodes;

			decrement_innerproduct(A, level, P, Q, permutation, p_G1, p_G2);
		}

		temp = permutation[level];
		permutation[level] = permutation[level + i];
		permutation[level + i] = temp;

	}

	return num_nodes;
}

int register_points_dfs(int num_points, double (*P)[3], double (*Q)[3], double max_rmsd, const treefilter_t* filter,
			uint8_t* best_permutation, double* p_rmsd, int* p_num_nodes_explored)
{
	double max_obj = max_rmsd * max_rmsd * num_points;

	uint8_t permutation[num_points];
	for (int i=0;i<num_points;i++)
		permutation[i] = i;

	double best_obj = INFINITY;
	double G1 = 0.0, G2 = 0.0, A[9] = {0};

	int num_nodes = recurse(num_points, 0, P, Q, max_obj, filter,
				permutation, &G1, &G2, A,
				best_permutation, &best_obj);
	*p_num_nodes_explored = num_nodes;

	if (best_obj < max_obj)
	{
		*p_rmsd = sqrt(best_obj / num_points);
		return 1;
	}
	else
	{
		*p_rmsd = INFINITY;
		return 0;
	}
}

}

