#include <cassert>
#include <cstring>
#include <cfloat>
#include <fstream>
#include <algorithm>
#include "reference_templates.hpp"
#include "convex_hull.hpp"
#include "canonical.hpp"
#include "graph_data.hpp"
#include "deformation_gradient.hpp"
#include "alloy_types.hpp"
#include "qcprot.hpp"
#include "quat.hpp"
#include "normalize_vertices.hpp"
#include "svdpolar/polar_decomposition.hpp"
#include "index_PTM.h"


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
	double normalized[15][3];
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
		max_degree = std::max(max_degree, degree[i]);

	return max_degree;
}

static void make_facets_clockwise(int num_facets, int8_t (*facets)[3], const double (*points)[3])
{
	double plane_normal[3];
	double origin[3] = {0, 0, 0};

	for (int i = 0;i<num_facets;i++)
		add_facet((double*)points, facets[i][0], facets[i][1], facets[i][2], facets[i], plane_normal, origin);
}

static void initialize_graphs(refdata_t* s)
{
	for (int i = 0;i<s->num_graphs;i++)
	{
		int8_t degree[s->num_nbrs];
		int _max_degree = graph_degree(s->num_facets, s->graphs[i].facets, s->num_nbrs, degree);
		assert(_max_degree <= s->max_degree);

		make_facets_clockwise(s->num_facets, s->graphs[i].facets, s->points);
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
				double scale,
				double (*normalized)[3],
				result_t* res)
{
	int num_points = s->num_nbrs + 1;
	const double (*ideal_points)[3] = s->points;
	int8_t inverse_labelling[num_points];
	int8_t forward_mapping[num_points];

	for (int i=0; i<num_points; i++)
		inverse_labelling[ canonical_labelling[i] ] = i;

	for (int i = 0;i<s->num_graphs;i++)
	{
		if (hash == s->graphs[i].hash)
		{
			graph_t* gref = &s->graphs[i];

			for (int j = 0;j<(int)gref->automorphisms.size();j++)
			{
				for (int k=0;k<num_points;k++)
					forward_mapping[gref->automorphisms[j][k]] = inverse_labelling[ gref->canonical_labelling[k] ];

				double A0[9], q[4], rmsd;
				double E0 = InnerProduct(A0, num_points, ideal_points, normalized, forward_mapping);
				FastCalcRMSDAndRotation(q, A0, &rmsd, E0, num_points, -1);

				if (rmsd < res->rmsd)
				{
					res->rmsd = rmsd;
					res->scale = scale;
					res->ref_struct = s;
					memcpy(res->q, q, 4 * sizeof(double));
					memcpy(res->mapping, forward_mapping, sizeof(int8_t) * num_points);
					memcpy(res->normalized, normalized, 3 * sizeof(double) * num_points);
				}
			}
		}
	}
}

static void match_general(refdata_t* s, double* points, result_t* res)
{
	double normalized[s->num_nbrs + 1][3];
	double scale = normalize_vertices(s->num_nbrs + 1, points, normalized);

	int8_t degree[s->num_nbrs];
	int8_t facets[s->num_facets][3];
	bool ok = get_convex_hull(s->num_nbrs + 1, (double*)normalized, s->num_facets, facets);
	if (!ok)
		return;

	int max_degree = graph_degree(s->num_facets, facets, s->num_nbrs, degree);
	if (max_degree > s->max_degree)
		return;

	if (s->type == PTM_MATCH_SC)
		for (int i = 0;i<s->num_nbrs;i++)
			if (degree[i] != 4)
				return;

	int8_t canonical_labelling[s->num_nbrs + 1];
	uint64_t hash = canonical_form(s->num_facets, facets, s->num_nbrs, degree, canonical_labelling);
	check_graphs(s, hash, canonical_labelling, scale, normalized, res);
}

static void match_fcc_hcp_ico(double* points, int32_t flags, result_t* res)
{
	int num_nbrs = structure_fcc.num_nbrs;
	int num_facets = structure_fcc.num_facets;
	int max_degree = structure_fcc.max_degree;

	double normalized[num_nbrs + 1][3];
	double scale = normalize_vertices(num_nbrs + 1, points, normalized);

	int8_t degree[num_nbrs];
	int8_t facets[num_facets][3];
	bool ok = get_convex_hull(num_nbrs + 1, (double*)normalized, num_facets, facets);
	if (!ok)
		return;

	int _max_degree = graph_degree(num_facets, facets, num_nbrs, degree);
	if (_max_degree > max_degree)
		return;

	int8_t canonical_labelling[num_nbrs + 1];
	uint64_t hash = canonical_form(num_facets, facets, num_nbrs, degree, canonical_labelling);
	if (flags & PTM_CHECK_FCC)	check_graphs(&structure_fcc, hash, canonical_labelling, scale, normalized, res);
	if (flags & PTM_CHECK_HCP)	check_graphs(&structure_hcp, hash, canonical_labelling, scale, normalized, res);
	if (flags & PTM_CHECK_ICO)	check_graphs(&structure_ico, hash, canonical_labelling, scale, normalized, res);
}

void index_PTM(	int num_points, double* points, int32_t* numbers, int32_t flags,
		int32_t* p_type, int32_t* p_alloy_type, double* p_scale, double* p_rmsd, double* q, double* F, double* F_res, double* U, double* P)
{
	if (flags & PTM_CHECK_SC)
		assert(num_points >= structure_sc.num_nbrs + 1);

	if (flags & PTM_CHECK_BCC)
		assert(num_points >= structure_bcc.num_nbrs + 1);

	if (flags & (PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO))
		assert(num_points >= structure_fcc.num_nbrs + 1);

	
	result_t res;
	res.ref_struct = NULL;
	res.rmsd = DBL_MAX;
	*p_type = PTM_MATCH_NONE;
	*p_alloy_type = PTM_ALLOY_NONE;

	if (flags & PTM_CHECK_SC)
		match_general(&structure_sc, points, &res);

	if (flags & (PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO))
		match_fcc_hcp_ico(points, flags, &res);

	if (flags & PTM_CHECK_BCC)
		match_general(&structure_bcc, points, &res);

	refdata_t* ref = res.ref_struct;
	if (ref != NULL)
	{
		*p_type = ref->type;

		if (numbers != NULL && ref->type == PTM_MATCH_FCC)
			*p_alloy_type = find_fcc_alloy_type(res.mapping, numbers);

		if (ref->type == PTM_MATCH_SC || ref->type == PTM_MATCH_FCC || ref->type == PTM_MATCH_BCC)
			rotate_quaternion_into_cubic_fundamental_zone(res.q);

		calculate_deformation_gradient(ref->num_nbrs + 1, ref->points, res.mapping, res.normalized, ref->penrose, F, F_res);
		left_sided_polar_decomposition_3x3(F, P, U);
	}

	*p_rmsd = res.rmsd;
	*p_scale = res.scale;
	memcpy(q, res.q, 4 * sizeof(double));
}

