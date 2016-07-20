#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cmath>
#include <cfloat>
#include <cassert>
#include "convex_hull_incremental.hpp"
#include "canonical.hpp"
#include "graph_data.hpp"
#include "deformation_gradient.hpp"
#include "alloy_types.hpp"
#include "neighbour_ordering.hpp"
#include "normalize_vertices.hpp"
#include "fundamental_mappings.hpp"
#include "qcprot/qcprot.hpp"
#include "qcprot/quat.hpp"
#include "polar_decomposition.hpp"
#include "index_ptm.h"


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
	const int8_t (*mapping)[15];
} refdata_t;

typedef struct
{
	double rmsd;
	double scale;
	double q[4];		//rotation in quaternion form (rigid body transformation)
	int8_t mapping[15];
	refdata_t* ref_struct;
} result_t;


//refdata_t structure_sc =  { .type = PTM_MATCH_SC,  .num_nbrs =  6, .num_facets =  8, .max_degree = 4, .num_graphs = NUM_SC_GRAPHS,  .graphs = graphs_sc,  .points = ptm_template_sc,  .penrose = penrose_sc , .mapping = mapping_sc };
refdata_t structure_sc =  { PTM_MATCH_SC,   6,  8, 4, NUM_SC_GRAPHS,  graphs_sc,  ptm_template_sc,  penrose_sc , mapping_sc };
refdata_t structure_fcc = { PTM_MATCH_FCC, 12, 20, 6, NUM_FCC_GRAPHS, graphs_fcc, ptm_template_fcc, penrose_fcc, mapping_fcc};
refdata_t structure_hcp = { PTM_MATCH_HCP, 12, 20, 6, NUM_HCP_GRAPHS, graphs_hcp, ptm_template_hcp, penrose_hcp, mapping_hcp};
refdata_t structure_ico = { PTM_MATCH_ICO, 12, 20, 6, NUM_ICO_GRAPHS, graphs_ico, ptm_template_ico, penrose_ico, mapping_ico};
refdata_t structure_bcc = { PTM_MATCH_BCC, 14, 24, 8, NUM_BCC_GRAPHS, graphs_bcc, ptm_template_bcc, penrose_bcc, mapping_bcc};

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

static int initialize_graphs(refdata_t* s)
{
	for (int i = 0;i<s->num_graphs;i++)
	{
		int8_t degree[PTM_MAX_NBRS];
		int _max_degree = graph_degree(s->num_facets, s->graphs[i].facets, s->num_nbrs, degree);
		assert(_max_degree <= s->max_degree);

		make_facets_clockwise(s->num_facets, s->graphs[i].facets, &s->points[1]);
		int ret = canonical_form(s->num_facets, s->graphs[i].facets, s->num_nbrs, degree, s->graphs[i].canonical_labelling, &s->graphs[i].hash);
		if (ret != 0)
			return ret;
	}

	return PTM_NO_ERROR;
}

bool ptm_initialized = false;
int ptm_initialize_global()
{
	if (ptm_initialized)
		return PTM_NO_ERROR;

	int ret = initialize_graphs(&structure_sc);
	ret |= initialize_graphs(&structure_fcc);
	ret |= initialize_graphs(&structure_hcp);
	ret |= initialize_graphs(&structure_ico);
	ret |= initialize_graphs(&structure_bcc);

	if (ret == PTM_NO_ERROR)
		ptm_initialized = true;

	return ret;
}

static void check_graphs(	refdata_t* s,
				uint64_t hash,
				int8_t* canonical_labelling,
				double (*normalized)[3],
				result_t* res)
{
	int num_points = s->num_nbrs + 1;
	const double (*ideal_points)[3] = s->points;
	int8_t inverse_labelling[PTM_MAX_POINTS];
	int8_t mapping[PTM_MAX_POINTS];

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
	int8_t degree[PTM_MAX_NBRS];
	int8_t facets[PTM_MAX_FACETS][3];
	int ret = get_convex_hull(s->num_nbrs + 1, (const double (*)[3])ch_points, s->num_facets, ch, facets);
	ch->ok = ret == 0;

#ifdef DEBUG
	printf("s->type: %d\tret: %d\n", s->type, ret);
#endif
	if (ret != 0)
		return PTM_NO_ERROR;

	int max_degree = graph_degree(s->num_facets, facets, s->num_nbrs, degree);
	if (max_degree > s->max_degree)
		return PTM_NO_ERROR;

	if (s->type == PTM_MATCH_SC)
		for (int i = 0;i<s->num_nbrs;i++)
			if (degree[i] != 4)
				return PTM_NO_ERROR;

	double normalized[PTM_MAX_POINTS][3];
	subtract_barycentre(s->num_nbrs + 1, points, normalized);

	int8_t canonical_labelling[PTM_MAX_POINTS];
	uint64_t hash = 0;
	ret = canonical_form(s->num_facets, facets, s->num_nbrs, degree, canonical_labelling, &hash);
	if (ret != PTM_NO_ERROR)
		return ret;

#ifdef DEBUG
	printf("hash: %lx\n", hash);
	printf("degree:\t");
	for (int i = 0;i<s->num_nbrs;i++)
		printf("%2d ", degree[i]);
	printf("\n");
#endif
	check_graphs(s, hash, canonical_labelling, normalized, res);
	return PTM_NO_ERROR;
}

static int match_fcc_hcp_ico(double (*ch_points)[3], double* points, int32_t flags, convexhull_t* ch, result_t* res)
{
	int num_nbrs = structure_fcc.num_nbrs;
	int num_facets = structure_fcc.num_facets;
	int max_degree = structure_fcc.max_degree;

	int8_t degree[PTM_MAX_NBRS];
	int8_t facets[PTM_MAX_FACETS][3];
	int ret = get_convex_hull(num_nbrs + 1, (const double (*)[3])ch_points, num_facets, ch, facets);
	ch->ok = ret == 0;

#ifdef DEBUG
	printf("s->type: %d\tret: %d\n", 2, ret);
#endif
	if (ret != 0)
		return PTM_NO_ERROR;

	int _max_degree = graph_degree(num_facets, facets, num_nbrs, degree);
	if (_max_degree > max_degree)
		return PTM_NO_ERROR;

	double normalized[PTM_MAX_POINTS][3];
	subtract_barycentre(num_nbrs + 1, points, normalized);

	int8_t canonical_labelling[PTM_MAX_POINTS];
	uint64_t hash = 0;
	ret = canonical_form(num_facets, facets, num_nbrs, degree, canonical_labelling, &hash);
	if (ret != PTM_NO_ERROR)
		return ret;

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
	return PTM_NO_ERROR;
}

/*static double calculate_lattice_constant(int type, double scale)
{
	assert(type >= 1 && type <= 5);
	double c[6] = {0, 2 / sqrt(2), 2 / sqrt(2), 14. / sqrt(3) - 7, 2 / sqrt(2), 1};
	return c[type] / scale;
}*/

static double calculate_lattice_constant(int type, double interatomic_distance)
{
	assert(type >= 1 && type <= 5);
	double c[6] = {0, 2 / sqrt(2), 2 / sqrt(2), 2. / sqrt(3), 2 / sqrt(2), 1};
	return c[type] * interatomic_distance;
}

static double calculate_interatomic_distance(int type, double scale)
{
	assert(type >= 1 && type <= 5);
	double c[6] = {0, 1, 1, (7. - 3.5 * sqrt(3)), 1, 1};
	return c[type] / scale;
}

int ptm_index(	ptm_local_handle_t local_handle, int num_points, double* unpermuted_points, int32_t* unpermuted_numbers, int32_t flags, bool topological_ordering,
		int32_t* p_type, int32_t* p_alloy_type, double* p_scale, double* p_rmsd, double* q, double* F, double* F_res, double* U, double* P, int8_t* mapping, double* p_interatomic_distance, double* p_lattice_constant)
{
	if (flags & PTM_CHECK_SC)
		assert(num_points >= structure_sc.num_nbrs + 1);

	if (flags & PTM_CHECK_BCC)
		assert(num_points >= structure_bcc.num_nbrs + 1);

	if (flags & (PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO))
		assert(num_points >= structure_fcc.num_nbrs + 1);


#define MAX_INPUT_POINTS 19
	assert(num_points <= MAX_INPUT_POINTS);

	int ret = 0;
	double ch_points[MAX_INPUT_POINTS][3];
	int8_t ordering[MAX_INPUT_POINTS];
	if (topological_ordering)
	{
		normalize_vertices(num_points, unpermuted_points, ch_points);
		ret = calculate_neighbour_ordering((void*)local_handle, num_points, (const double (*)[3])ch_points, ordering);
		if (ret != 0)
			topological_ordering = false;
	}

	if (!topological_ordering)
		for (int i=0;i<num_points;i++)
			ordering[i] = i;

	double points[PTM_MAX_POINTS][3];
	int32_t numbers[PTM_MAX_POINTS];
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
	if (p_alloy_type != NULL)
		*p_alloy_type = PTM_ALLOY_NONE;

	if (mapping != NULL)
		memset(mapping, -1, num_points * sizeof(int8_t));


	if (flags & PTM_CHECK_SC)
	{
		ret = match_general(&structure_sc, ch_points, (double*)points, &ch, &res);
		//if (ret != PTM_NO_ERROR)
		//	return ret;
#ifdef DEBUG
		printf("match sc  ret: %d\t%p\n", ret, res.ref_struct);
#endif
	}

	if (flags & (PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO))
	{
		ret = match_fcc_hcp_ico(ch_points, (double*)points, flags, &ch, &res);
		//if (ret != PTM_NO_ERROR)
		//	return ret;
#ifdef DEBUG
		printf("match fcc ret: %d\t%p\n", ret, res.ref_struct);
#endif
	}

	if (flags & PTM_CHECK_BCC)
	{
		ret = match_general(&structure_bcc, ch_points, (double*)points, &ch, &res);
		//if (ret != PTM_NO_ERROR)
		//	return ret;
#ifdef DEBUG
		printf("match bcc ret: %d\t%p\n", ret, res.ref_struct);
#endif
	}

	refdata_t* ref = res.ref_struct;
	if (ref != NULL)
	{
		*p_type = ref->type;

		if (p_alloy_type != NULL && unpermuted_numbers != NULL)
		{
			if (ref->type == PTM_MATCH_FCC)
				*p_alloy_type = find_fcc_alloy_type(res.mapping, numbers);
			else if (ref->type == PTM_MATCH_BCC)
				*p_alloy_type = find_bcc_alloy_type(res.mapping, numbers);
		}

		int bi = -1;
		if      (ref->type == PTM_MATCH_SC)	bi = rotate_quaternion_into_cubic_fundamental_zone(res.q);
		else if (ref->type == PTM_MATCH_FCC)	bi = rotate_quaternion_into_cubic_fundamental_zone(res.q);
		else if (ref->type == PTM_MATCH_BCC)	bi = rotate_quaternion_into_cubic_fundamental_zone(res.q);
		else if (ref->type == PTM_MATCH_ICO)	bi = rotate_quaternion_into_icosahedral_fundamental_zone(res.q);
		else if (ref->type == PTM_MATCH_HCP)	bi = rotate_quaternion_into_hcp_fundamental_zone(res.q);

		int8_t temp[15];
		for (int i=0;i<ref->num_nbrs+1;i++)
			temp[ref->mapping[bi][i]] = res.mapping[i];

		memcpy(res.mapping, temp, (ref->num_nbrs+1) * sizeof(int8_t));

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
				polar_decomposition_3x3(F, false, U, P);
		}

		if (mapping != NULL)
			for (int i=0;i<ref->num_nbrs + 1;i++)
				mapping[i] = ordering[res.mapping[i]];

		double interatomic_distance = calculate_interatomic_distance(ref->type, res.scale);
		double lattice_constant = calculate_lattice_constant(ref->type, interatomic_distance);

		if (p_interatomic_distance != NULL)
			*p_interatomic_distance = interatomic_distance;

		if (p_lattice_constant != NULL)
			*p_lattice_constant = lattice_constant;
	}

	*p_rmsd = res.rmsd;
	*p_scale = res.scale;
	memcpy(q, res.q, 4 * sizeof(double));

	return PTM_NO_ERROR;
}

ptm_local_handle_t ptm_initialize_local()
{
	return (ptm_local_handle_t)voronoi_initialize_local();
}

void ptm_uninitialize_local(ptm_local_handle_t ptr)
{
	voronoi_uninitialize_local(ptr);
}

