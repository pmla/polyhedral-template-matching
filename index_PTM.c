#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include "reference_templates.h"
#include "convex_hull_incremental.h"
#include "canonical.h"
#include "graph_data.h"
#include "deformation_gradient.h"
#include "alloy_types.h"
#include "qcprot.h"
#include "quat.h"
#include "normalize_vertices.h"
#include "svdpolar/polar_decomposition.h"
#include "neighbour_ordering.h"
#include "index_PTM.h"


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))


typedef struct
{
	int type;
	int num_nbrs;
	int num_facets;
	int max_degree;
	int num_graphs;
	graph_t* graphs;
	const double (*points)[3];
	const double (*penrose)[3];
} refdata_t;

typedef struct
{
	double rmsd;
	double scale;
	double q[4];		//rotation in quaternion form (rigid body transformation)
	int8_t mapping[15];
	refdata_t* ref_struct;
} result_t;

refdata_t structure_sc =  { .type = PTM_MATCH_SC,  .num_nbrs =  6, .num_facets =  8, .max_degree = 4, .num_graphs = NUM_SC_GRAPHS,  .graphs = graphs_sc,  .points = points_sc,  .penrose = penrose_sc };
refdata_t structure_fcc = { .type = PTM_MATCH_FCC, .num_nbrs = 12, .num_facets = 20, .max_degree = 6, .num_graphs = NUM_FCC_GRAPHS, .graphs = graphs_fcc, .points = points_fcc, .penrose = penrose_fcc};
refdata_t structure_hcp = { .type = PTM_MATCH_HCP, .num_nbrs = 12, .num_facets = 20, .max_degree = 6, .num_graphs = NUM_HCP_GRAPHS, .graphs = graphs_hcp, .points = points_hcp, .penrose = penrose_hcp};
refdata_t structure_ico = { .type = PTM_MATCH_ICO, .num_nbrs = 12, .num_facets = 20, .max_degree = 6, .num_graphs = NUM_ICO_GRAPHS, .graphs = graphs_ico, .points = points_ico, .penrose = penrose_ico};
refdata_t structure_bcc = { .type = PTM_MATCH_BCC, .num_nbrs = 14, .num_facets = 24, .max_degree = 8, .num_graphs = NUM_BCC_GRAPHS, .graphs = graphs_bcc, .points = points_bcc, .penrose = penrose_bcc};

static int graph_degree(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree)
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
		max_degree = MAX(max_degree, degree[i]);

	return max_degree;
}

static void make_facets_clockwise(int num_facets, int8_t (*facets)[3], const double (*points)[3])
{
	double plane_normal[3];
	double origin[3] = {0, 0, 0};

	for (int i = 0;i<num_facets;i++)
		add_facet(points, facets[i][0], facets[i][1], facets[i][2], facets[i], plane_normal, origin);
}

static void initialize_graphs(refdata_t* s)
{
	for (int i = 0;i<s->num_graphs;i++)
	{
		int8_t degree[s->num_nbrs];
		int _max_degree = graph_degree(s->num_facets, s->graphs[i].facets, s->num_nbrs, degree);
		assert(_max_degree <= s->max_degree);

		make_facets_clockwise(s->num_facets, s->graphs[i].facets, &s->points[1]);
		s->graphs[i].hash = canonical_form(s->num_facets, s->graphs[i].facets, s->num_nbrs, degree, s->graphs[i].canonical_labelling);
	}
}

int initialize_PTM()
{
	initialize_graphs(&structure_sc);
	initialize_graphs(&structure_fcc);
	initialize_graphs(&structure_hcp);
	initialize_graphs(&structure_ico);
	initialize_graphs(&structure_bcc);
	return 0;
}

static void check_graphs(	refdata_t* s,
				uint64_t hash,
				int8_t* canonical_labelling,
				double (*normalized)[3],
				result_t* res)
{
	int num_points = s->num_nbrs + 1;
	const double (*ideal_points)[3] = s->points;
	int8_t inverse_labelling[num_points];
	int8_t mapping[num_points];

	for (int i=0; i<num_points; i++)
		inverse_labelling[ canonical_labelling[i] ] = i;

	double G1 = 0, G2 = 0;
	for (int i=0;i<num_points;i++)
	{
		double x1 = ideal_points[i][0];
		double y1 = ideal_points[i][1];
		double z1 = ideal_points[i][2];

		double x2 = normalized[i][0];
		double y2 = normalized[i][1];
		double z2 = normalized[i][2];

		G1 += x1 * x1 + y1 * y1 + z1 * z1;
		G2 += x2 * x2 + y2 * y2 + z2 * z2;
	}
	double E0 = (G1 + G2) / 2;

	for (int i = 0;i<s->num_graphs;i++)
	{
//printf("graph hash: %lx\n", s->graphs[i].hash);
		if (hash == s->graphs[i].hash)
		{
			graph_t* gref = &s->graphs[i];

			for (int j = 0;j<gref->num_automorphisms;j++)
			{
				for (int k=0;k<num_points;k++)
					mapping[automorphisms[gref->automorphism_index + j][k]] = inverse_labelling[ gref->canonical_labelling[k] ];

				double A0[9], q[4], rmsd;
				InnerProduct(A0, num_points, ideal_points, normalized, mapping);

				double rot[9];
				FastCalcRMSDAndRotation(q, A0, &rmsd, E0, num_points, -1, rot);

				double k0 = 0;
				for (int ii=0;ii<num_points;ii++)
				{
					for (int jj=0;jj<3;jj++)
					{
						double v = 0.0;
						for (int kk=0;kk<3;kk++)
							v += rot[jj*3+kk] * ideal_points[ii][kk];

						k0 += v * normalized[mapping[ii]][jj];
					}
				}

				double scale = k0 / G2;
				rmsd = sqrt(fabs(G1 - scale*k0) / num_points);
//printf("got one: %f\n", rmsd);

				if (rmsd < res->rmsd)
				{
					res->rmsd = rmsd;
					res->scale = scale;
					res->ref_struct = s;
					memcpy(res->q, q, 4 * sizeof(double));
					memcpy(res->mapping, mapping, sizeof(int8_t) * num_points);
				}
			}
		}
	}
}

static int match_general(refdata_t* s, double (*ch_points)[3], double* points, convexhull_t* ch, result_t* res)
{
	int8_t degree[s->num_nbrs];
	int8_t facets[s->num_facets][3];
	int ret = get_convex_hull(s->num_nbrs + 1, (const double (*)[3])ch_points, s->num_facets, ch, facets);
	ch->ok = ret == 0;

#ifdef DEBUG
	printf("s->type: %d\tret: %d\n", s->type, ret);
#endif
	if (ret != 0)
		return ret;

	int max_degree = graph_degree(s->num_facets, facets, s->num_nbrs, degree);
	if (max_degree > s->max_degree)
		return -7;

	if (s->type == PTM_MATCH_SC)
		for (int i = 0;i<s->num_nbrs;i++)
			if (degree[i] != 4)
				return -8;

	double normalized[s->num_nbrs + 1][3];
	subtract_barycentre(s->num_nbrs + 1, points, normalized);

	int8_t canonical_labelling[s->num_nbrs + 1];
	uint64_t hash = canonical_form(s->num_facets, facets, s->num_nbrs, degree, canonical_labelling);
#ifdef DEBUG
	printf("hash: %lx\n", hash);
	printf("degree:\t");
	for (int i = 0;i<s->num_nbrs;i++)
		printf("%2d ", degree[i]);
	printf("\n");
#endif
	check_graphs(s, hash, canonical_labelling, normalized, res);
	return 0;
}

static int match_fcc_hcp_ico(double (*ch_points)[3], double* points, int32_t flags, convexhull_t* ch, result_t* res)
{
	int num_nbrs = structure_fcc.num_nbrs;
	int num_facets = structure_fcc.num_facets;
	int max_degree = structure_fcc.max_degree;

	int8_t degree[num_nbrs];
	int8_t facets[num_facets][3];
	int ret = get_convex_hull(num_nbrs + 1, (const double (*)[3])ch_points, num_facets, ch, facets);
	ch->ok = ret == 0;

#ifdef DEBUG
	printf("s->type: %d\tret: %d\n", 2, ret);
#endif
	if (ret != 0)
		return ret;

	int _max_degree = graph_degree(num_facets, facets, num_nbrs, degree);
	if (_max_degree > max_degree)
		return -9;

	double normalized[num_nbrs + 1][3];
	subtract_barycentre(num_nbrs + 1, points, normalized);

	int8_t canonical_labelling[num_nbrs + 1];
	uint64_t hash = canonical_form(num_facets, facets, num_nbrs, degree, canonical_labelling);
#ifdef DEBUG
	printf("hash: %lx\n", hash);
	printf("degree:\t");
	for (int i = 0;i<num_nbrs;i++)
		printf("%2d ", degree[i]);
	printf("\n");
#endif
	if (flags & PTM_CHECK_FCC)	check_graphs(&structure_fcc, hash, canonical_labelling, normalized, res);
	if (flags & PTM_CHECK_HCP)	check_graphs(&structure_hcp, hash, canonical_labelling, normalized, res);
	if (flags & PTM_CHECK_ICO)	check_graphs(&structure_ico, hash, canonical_labelling, normalized, res);
	return 0;
}

#ifdef DEBUG
int failcount = 0;
#endif

void index_PTM(	int num_points, double* unpermuted_points, int32_t* unpermuted_numbers, int32_t flags, bool topological_ordering,
		int32_t* p_type, int32_t* p_alloy_type, double* p_scale, double* p_rmsd, double* q, double* F, double* F_res, double* U, double* P)
{
	if (flags & PTM_CHECK_SC)
		assert(num_points >= structure_sc.num_nbrs + 1);

	if (flags & PTM_CHECK_BCC)
		assert(num_points >= structure_bcc.num_nbrs + 1);

	if (flags & (PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO))
		assert(num_points >= structure_fcc.num_nbrs + 1);


	assert(num_points <= 19);

	int ret = 0;
	double ch_points[19][3];
	int8_t ordering[19];
	if (topological_ordering)
	{
		normalize_vertices(num_points, unpermuted_points, ch_points);
		ret = calculate_neighbour_ordering(num_points, (const double (*)[3])ch_points, ordering);
		if (ret != 0)
		{
#ifdef DEBUG
			failcount++;
#endif
			topological_ordering = false;
		}
	}

	if (!topological_ordering)
		for (int i=0;i<num_points;i++)
			ordering[i] = i;

	double points[15][3];
	int32_t numbers[19];
	num_points = MIN(15, num_points);
	for (int i=0;i<num_points;i++)
	{
		memcpy(points[i], &unpermuted_points[3 * ordering[i]], 3 * sizeof(double));

		if (unpermuted_numbers != NULL)
			numbers[i] = unpermuted_numbers[ordering[i]];
	}

	convexhull_t ch;
	ch.ok = false;
	normalize_vertices(num_points, (double*)points, ch_points);

#ifdef DEBUG
	for (int i = 0;i<num_points;i++)
		printf("%.2f\t%.2f\t%.2f\n", unpermuted_points[i*3 + 0], unpermuted_points[i*3 + 1], unpermuted_points[i*3 + 2]);
#endif

	result_t res;
	res.ref_struct = NULL;
	res.rmsd = INFINITY;
	*p_type = PTM_MATCH_NONE;
	*p_alloy_type = PTM_ALLOY_NONE;

	if (flags & PTM_CHECK_SC)
	{
		ret = match_general(&structure_sc, ch_points, (double*)points, &ch, &res);
#ifdef DEBUG
		printf("match sc  ret: %d\t%p\n", ret, res.ref_struct);
#endif
	}

	if (flags & (PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO))
	{
		ret = match_fcc_hcp_ico(ch_points, (double*)points, flags, &ch, &res);
#ifdef DEBUG
		printf("match fcc ret: %d\t%p\n", ret, res.ref_struct);
#endif
	}

	if (flags & PTM_CHECK_BCC)
	{
		ret = match_general(&structure_bcc, ch_points, (double*)points, &ch, &res);
#ifdef DEBUG
		printf("match bcc ret: %d\t%p\n", ret, res.ref_struct);
#endif
	}

	refdata_t* ref = res.ref_struct;
	if (ref != NULL)
	{
		*p_type = ref->type;

		if (unpermuted_numbers != NULL)
		{
			if (ref->type == PTM_MATCH_FCC)
				*p_alloy_type = find_fcc_alloy_type(res.mapping, numbers);
			else if (ref->type == PTM_MATCH_BCC)
				*p_alloy_type = find_bcc_alloy_type(res.mapping, numbers);
		}

		if (ref->type == PTM_MATCH_SC || ref->type == PTM_MATCH_FCC || ref->type == PTM_MATCH_BCC)
			rotate_quaternion_into_cubic_fundamental_zone(res.q);

		if (F != NULL && F_res != NULL)
		{
			subtract_barycentre(ref->num_nbrs + 1, (double*)points, ch_points);
			for (int i = 0;i<ref->num_nbrs + 1;i++)
			{
				ch_points[i][0] *= res.scale;
				ch_points[i][1] *= res.scale;
				ch_points[i][2] *= res.scale;
			}
			calculate_deformation_gradient(ref->num_nbrs + 1, ref->points, res.mapping, ch_points, ref->penrose, F, F_res);

			if (P != NULL && U != NULL)
				left_sided_polar_decomposition_3x3(F, P, U);
		}
	}

	*p_rmsd = res.rmsd;
	*p_scale = res.scale;
	memcpy(q, res.q, 4 * sizeof(double));
}

