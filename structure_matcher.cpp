#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <algorithm>
#include "convex_hull_incremental.hpp"
#include "canonical.hpp"
#include "graph_data.hpp"
#include "normalize_vertices.hpp"
#include "qcprot/polar.hpp"
#include "initialize_data.hpp"
#include "structure_matcher.hpp"
#include "ptm_constants.h"


extern refdata_t structure_sc;
extern refdata_t structure_fcc;
extern refdata_t structure_hcp;
extern refdata_t structure_ico;
extern refdata_t structure_bcc;


static double calc_rmsd(int num_points, const double (*ideal_points)[3], double (*normalized)[3], int8_t* mapping,
			double G1, double G2, double E0, double* q, double* p_scale)
{
	double A0[9];
	InnerProduct(A0, num_points, ideal_points, normalized, mapping);

	double nrmsdsq, rot[9];
	FastCalcRMSDAndRotation(A0, E0, &nrmsdsq, q, rot);

	double k0 = 0;
	for (int i=0;i<num_points;i++)
	{
		for (int j=0;j<3;j++)
		{
			double v = 0.0;
			for (int k=0;k<3;k++)
				v += rot[j*3+k] * ideal_points[i][k];

			k0 += v * normalized[mapping[i]][j];
		}
	}

	double scale = k0 / G2;
	*p_scale = scale;
	return sqrt(fabs(G1 - scale*k0) / num_points);
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
		if (hash != s->graphs[i].hash)
			continue;

		graph_t* gref = &s->graphs[i];
		for (int j = 0;j<gref->num_automorphisms;j++)
		{
			for (int k=0;k<num_points;k++)
				mapping[automorphisms[gref->automorphism_index + j][k]] = inverse_labelling[ gref->canonical_labelling[k] ];

			double q[4], scale = 0;
			double rmsd = calc_rmsd(num_points, ideal_points, normalized, mapping, G1, G2, E0, q, &scale);
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

int match_general(refdata_t* s, double (*ch_points)[3], double (*points)[3], convexhull_t* ch, result_t* res)
{
	int8_t degree[PTM_MAX_NBRS];
	int8_t facets[PTM_MAX_FACETS][3];

	int ret = get_convex_hull(s->num_nbrs + 1, (const double (*)[3])ch_points, s->num_facets, ch, facets);
	ch->ok = ret == 0;
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

	check_graphs(s, hash, canonical_labelling, normalized, res);
	return PTM_NO_ERROR;
}

int match_fcc_hcp_ico(double (*ch_points)[3], double (*points)[3], int32_t flags, convexhull_t* ch, result_t* res)
{
	int num_nbrs = structure_fcc.num_nbrs;
	int num_facets = structure_fcc.num_facets;
	int max_degree = structure_fcc.max_degree;

	int8_t degree[PTM_MAX_NBRS];
	int8_t facets[PTM_MAX_FACETS][3];

	int ret = get_convex_hull(num_nbrs + 1, (const double (*)[3])ch_points, num_facets, ch, facets);
	ch->ok = ret == 0;
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

	if (flags & PTM_CHECK_FCC)	check_graphs(&structure_fcc, hash, canonical_labelling, normalized, res);
	if (flags & PTM_CHECK_HCP)	check_graphs(&structure_hcp, hash, canonical_labelling, normalized, res);
	if (flags & PTM_CHECK_ICO)	check_graphs(&structure_ico, hash, canonical_labelling, normalized, res);
	return PTM_NO_ERROR;
}

