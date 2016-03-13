#include <cassert>
#include <cstring>
#include <cfloat>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "reference_templates.hpp"
#include "convex_hull.hpp"
#include "canonical.hpp"
#include "graph_data.hpp"
#include "index_PTM.h"


static double frand(double fmin, double fmax)
{
	double f = (double)rand() / RAND_MAX;
	return fmin + f * (fmax - fmin);
}

static void random_permutation(int n, int* permutation)
{
	int range[n];
	for (int i = 0;i<n;i++)
		range[i] = i;

	int num_left = n;
	for (int i = 0;i<n;i++)
	{
		int index = rand() % num_left;
		permutation[i] = range[index];
		range[index] = range[--num_left];
	}
}


static void random_quaternion(double* q)
{
	double x0 = frand(0.0, 1.0);
	double x1 = frand(0.0, 2 * M_PI);
	double x2 = frand(0.0, 2 * M_PI);

	double r1 = sqrt(1 - x0);
	double r2 = sqrt(x0);

	double s1 = sin(x1);
	double c1 = cos(x1);
	double s2 = sin(x2);
	double c2 = cos(x2);

	q[0] = s1*r1;
	q[1] = c1*r1;
	q[2] = s2*r2;
	q[3] = c2*r2;
}

static void quaternion_to_rotation_matrix(double* q, double* u)
{
	double a = q[0];
	double b = q[1];
	double c = q[2];
	double d = q[3];

	u[0] = a*a + b*b - c*c - d*d;
	u[1] = 2*b*c - 2*a*d;
	u[2] = 2*b*d + 2*a*c;

	u[3] = 2*b*c + 2*a*d;
	u[4] = a*a - b*b + c*c - d*d;
	u[5] = 2*c*d - 2*a*b;

	u[6] = 2*b*d - 2*a*c;
	u[7] = 2*c*d + 2*a*b;
	u[8] = a*a - b*b - c*c + d*d;
}

static void random_rotation_matrix(double* U)
{
	double q[4];
	random_quaternion(q);
	quaternion_to_rotation_matrix(q, U);
}

static void rotate_vector(const double* U, const double* x, double* b)
{
	b[0] = U[0] * x[0] + U[1] * x[1] + U[2] * x[2];
	b[1] = U[3] * x[0] + U[4] * x[1] + U[5] * x[2];
	b[2] = U[6] * x[0] + U[7] * x[1] + U[8] * x[2];
}

typedef struct
{
	int type;
	int num_nbrs;
	int num_facets;
	int max_degree;
	int num_graphs;
	graph_t* graphs;
	const double (*points)[3];
} _refdata_t;

_refdata_t _structure_fcc = { .type = PTM_MATCH_FCC, .num_nbrs = 12, .num_facets = 20, .max_degree = 6, .num_graphs = NUM_FCC_GRAPHS, .graphs = graphs_fcc, .points = points_fcc};
_refdata_t _structure_hcp = { .type = PTM_MATCH_HCP, .num_nbrs = 12, .num_facets = 20, .max_degree = 6, .num_graphs = NUM_HCP_GRAPHS, .graphs = graphs_hcp, .points = points_hcp};

int self_test()
{
	int num_nbrs = _structure_fcc.num_nbrs;
	int num_facets = _structure_fcc.num_facets;
	int max_degree = _structure_fcc.max_degree;

	int permutation[num_nbrs+1];
	random_permutation(num_nbrs, permutation);
	permutation[num_nbrs] = num_nbrs;

	double U[9];
	random_rotation_matrix(U);

	int num_points = num_nbrs + 1;
	double points[num_points][3];
	for (int i = 0;i<num_points;i++)
		rotate_vector(U, _structure_hcp.points[i], points[permutation[i]]);

	int8_t facets[num_facets][3];
	bool ok = get_convex_hull(num_nbrs + 1, (double*)points, num_facets, facets);
	if (!ok)
	{
		for (int i = 0;i<num_points;i++)
			printf("%-3d ", permutation[i]);
		printf("\n");
	}
	return ok;
	//if (!ok)
	//	return 0;

	/*int8_t degree[num_nbrs];
	int _max_degree = graph_degree(num_facets, facets, num_nbrs, degree);
	if (_max_degree > max_degree)
		return 0;

	make_facets_clockwise(num_facets, facets, normalized);

	int8_t canonical_labelling[num_nbrs + 1];
	uint64_t hash = canonical_form(num_facets, facets, num_nbrs, max_degree, degree, canonical_labelling);
	if (flags & PTM_CHECK_FCC)	num_matches += check_graphs(NULL, &structure_fcc, hash, canonical_labelling, normalized, numbers, &res[num_matches]) ? 1 : 0;
	if (flags & PTM_CHECK_HCP)	num_matches += check_graphs(NULL, &structure_hcp, hash, canonical_labelling, normalized,    NULL, &res[num_matches]) ? 1 : 0;
	if (flags & PTM_CHECK_ICO)	num_matches += check_graphs(NULL, &structure_ico, hash, canonical_labelling, normalized,    NULL, &res[num_matches]) ? 1 : 0;
	return 0;*/
}

