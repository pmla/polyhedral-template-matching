/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "ptm_alloy_types.h"
#include "ptm_constants.h"
#include "ptm_convex_hull_incremental.h"
#include "ptm_deformation_gradient.h"
#include "ptm_functions.h"
#include "ptm_graph_data.h"
#include "ptm_initialize_data.h"
#include "ptm_neighbour_ordering.h"
#include "ptm_normalize_vertices.h"
#include "ptm_polar.h"
#include "ptm_quat.h"
#include "ptm_structure_matcher.h"
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>


static double calculate_interatomic_distance(int type, double scale) {
	assert(type >= 1 && type <= 8);

	// these values should be equal to norm(template[1])
	double c[9] = {0,
		       1,
		       1,
		       (7. - 3.5 * sqrt(3)),
		       1,
		       1,
		       sqrt(3) * 4. / (6 * sqrt(2) + sqrt(3)),
		       sqrt(3) * 4. / (6 * sqrt(2) + sqrt(3)),
		       -3. / 11 + 6 * sqrt(3) / 11};
	return c[type] / scale;
}

static double calculate_lattice_constant(int type,
					 double interatomic_distance) {
	assert(type >= 1 && type <= 8);
	double c[9] = {0, 2 / sqrt(2), 2 / sqrt(2), 2. / sqrt(3), 2 / sqrt(2),
		       1, 4 / sqrt(3), 4 / sqrt(3), sqrt(3)};
	return c[type] * interatomic_distance;
}

static int rotate_into_fundamental_zone(int type,
					bool output_conventional_orientation,
					double *q) {
	if (type == PTM_MATCH_SC)
		return ptm::rotate_quaternion_into_cubic_fundamental_zone(q);
	if (type == PTM_MATCH_FCC)
		return ptm::rotate_quaternion_into_cubic_fundamental_zone(q);
	if (type == PTM_MATCH_BCC)
		return ptm::rotate_quaternion_into_cubic_fundamental_zone(q);
	if (type == PTM_MATCH_ICO)
		return ptm::rotate_quaternion_into_icosahedral_fundamental_zone(q);

	if (type == PTM_MATCH_HCP || type == PTM_MATCH_GRAPHENE) {
		if (output_conventional_orientation) {
			return ptm::rotate_quaternion_into_hcp_conventional_fundamental_zone(q);
		} else {
			return ptm::rotate_quaternion_into_hcp_fundamental_zone(q);
		}
	}

	if (type == PTM_MATCH_DCUB) {
		if (output_conventional_orientation) {
			return ptm::rotate_quaternion_into_cubic_fundamental_zone(q);
		} else {
			return ptm::rotate_quaternion_into_diamond_cubic_fundamental_zone(q);
		}
	}

	if (type == PTM_MATCH_DHEX) {
		if (output_conventional_orientation) {
			return ptm::rotate_quaternion_into_hcp_conventional_fundamental_zone(q);
		} else {
			return ptm::rotate_quaternion_into_diamond_hexagonal_fundamental_zone(q);
		}
	}

	return -1;
}

int ptm_remap_template(	int type, bool output_conventional_orientation, double* qtarget, double* q,
			double* qmapped, double* p_disorientation, int8_t* mapping)
{
	if (type == PTM_MATCH_NONE)
		return -1;

	const ptm::refdata_t* refs[] = {	NULL,
						&ptm::structure_fcc,
						&ptm::structure_hcp,
						&ptm::structure_bcc,
						&ptm::structure_ico,
						&ptm::structure_sc,
						&ptm::structure_dcub,
						&ptm::structure_dhex,
						&ptm::structure_graphene	};

	const ptm::refdata_t* ref = refs[type];

	int8_t temp[PTM_MAX_POINTS];
	memset(temp, -1, PTM_MAX_POINTS * sizeof(int8_t));

	double qrel[4];
	if (qtarget != NULL) {
		double qinv[4] = {-qtarget[0], qtarget[1], qtarget[2], qtarget[3]};
		ptm::quat_rot(q, qinv, qrel);
	}
	else
	{
		memcpy(qrel, q, 4 * sizeof(double));
	}

	int bi = rotate_into_fundamental_zone(ref->type, output_conventional_orientation, qrel);
	if (bi < 0)
		return bi;

	if (output_conventional_orientation & (    ref->type == PTM_MATCH_HCP
						|| ref->type == PTM_MATCH_GRAPHENE
						|| ref->type == PTM_MATCH_DCUB
						|| ref->type == PTM_MATCH_DHEX))
	{
		for (int i=0;i<ref->num_nbrs+1;i++)
			temp[ref->mapping_conventional[bi][i]] = mapping[i];
	}
	else
	{
		for (int i=0;i<ref->num_nbrs+1;i++)
			temp[ref->mapping[bi][i]] = mapping[i];
	}

	if (qtarget != NULL) {
		double identity[4] = {1, 0, 0, 0};
		*p_disorientation = ptm::quat_misorientation(identity, qrel);
		ptm::quat_rot(qrel, q, qmapped);
	}
	else {
		memcpy(q, qrel, 4 * sizeof(double));
	}

	memcpy(mapping, temp, (ref->num_nbrs + 1) * sizeof(int8_t));
	return bi;
}

static void output_data(ptm::result_t *res, ptm::atomicenv_t* env,
			bool output_conventional_orientation, int32_t *p_type,
			int32_t *p_alloy_type, double *p_scale, double *p_rmsd,
			double *q, double *F, double *F_res, double *U,
			double *P, double *p_interatomic_distance,
			double *p_lattice_constant,
			int* p_best_template_index, const double (**p_best_template)[3],
			int8_t *output_indices) {
	const ptm::refdata_t *ref = res->ref_struct;
	if (ref == NULL)
		return;

	*p_type = ref->type;
	if (p_alloy_type != NULL)
		*p_alloy_type = ptm::find_alloy_type(ref, res->mapping, env->numbers);

	int bi = ptm_remap_template(	ref->type, output_conventional_orientation,
					NULL, res->q, NULL, NULL, res->mapping);
	if (bi < 0)
		return;

	int best_template_index = 0;
	const double (*ref_template)[3] = ref->points;
	const double (*ref_penrose)[3] = ref->penrose;
	if (output_conventional_orientation & (    ref->type == PTM_MATCH_HCP
						|| ref->type == PTM_MATCH_GRAPHENE
						|| ref->type == PTM_MATCH_DCUB
						|| ref->type == PTM_MATCH_DHEX))
	{
		if (ref->template_indices[bi] == 1)
		{
			ref_template = ref->points_alt1;
			ref_penrose = ref->penrose_alt1;
		}
		else if (ref->template_indices[bi] == 2)
		{
			ref_template = ref->points_alt2;
			ref_penrose = ref->penrose_alt2;
		}
		else if (ref->template_indices[bi] == 3)
		{
			ref_template = ref->points_alt3;
			ref_penrose = ref->penrose_alt3;
		}

		best_template_index = ref->template_indices[bi];
	}

	if (p_best_template_index != NULL)
		*p_best_template_index = best_template_index;

	if (p_best_template != NULL)
		*p_best_template = ref_template;

	if (F != NULL && F_res != NULL) {
		double scaled_points[PTM_MAX_INPUT_POINTS][3];

		ptm::subtract_barycentre(ref->num_nbrs + 1, env->points, scaled_points);
		for (int i = 0; i < ref->num_nbrs + 1; i++) {
			scaled_points[i][0] *= res->scale;
			scaled_points[i][1] *= res->scale;
			scaled_points[i][2] *= res->scale;
		}

		ptm::calculate_deformation_gradient(ref->num_nbrs + 1, ref_template,
						    res->mapping, scaled_points, ref_penrose,
						    F, F_res);
		if (ref->type == PTM_MATCH_GRAPHENE) // hack for pseudo-2d structures
			F[8] = 1;

		if (P != NULL && U != NULL)
			ptm::polar_decomposition_3x3(F, false, U, P);
	}

	if (output_indices != NULL)
		for (int i = 0; i < ref->num_nbrs + 1; i++)
			output_indices[i] = env->ordering[res->mapping[i]];

	double interatomic_distance = calculate_interatomic_distance(ref->type, res->scale);
	double lattice_constant = calculate_lattice_constant(ref->type, interatomic_distance);

	if (p_interatomic_distance != NULL)
		*p_interatomic_distance = interatomic_distance;

	if (p_lattice_constant != NULL)
		*p_lattice_constant = lattice_constant;

	*p_rmsd = res->rmsd;
	*p_scale = res->scale;
	memcpy(q, res->q, 4 * sizeof(double));
}

extern bool ptm_initialized;

int ptm_index(ptm_local_handle_t local_handle,
	      int num_input_points, double (*input_points)[3], int32_t* input_numbers,
	      int32_t flags,
	      bool output_conventional_orientation, int32_t *p_type,
	      int32_t *p_alloy_type, double *p_scale, double *p_rmsd, double *q,
	      double *F, double *F_res, double *U, double *P,
	      double *p_interatomic_distance, double *p_lattice_constant,
	      int* p_best_template_index, const double (**p_best_template)[3],
	      int8_t *output_indices)
{
	int ret = 0;
	assert(ptm_initialized);
	if (!ptm_initialized)	//assert is not active in OVITO release build
		return -1;

	assert(num_input_points <= PTM_MAX_INPUT_POINTS);

	//-------- initialize output values with failure case --------
	ptm::result_t res;
	res.ref_struct = NULL;
	res.rmsd = INFINITY;

	if (output_indices != NULL)
		memset(output_indices, -1, PTM_MAX_INPUT_POINTS * sizeof(int8_t));

	*p_type = PTM_MATCH_NONE;
	if (p_alloy_type != NULL)
		*p_alloy_type = PTM_ALLOY_NONE;
	//------------------------------------------------------------

	int32_t dummy_input_numbers[PTM_MAX_INPUT_POINTS] = {0};
	if (input_numbers == NULL) {
		input_numbers = dummy_input_numbers;
		if (p_alloy_type != NULL)
			*p_alloy_type = PTM_ALLOY_NONE;
		p_alloy_type = NULL;
	}

	ptm::atomicenv_t env, dmn_env, grp_env;

	ptm::convexhull_t ch;
	double ch_points[PTM_MAX_INPUT_POINTS][3];

	ptm::calculate_neighbour_ordering(local_handle, num_input_points, input_points, input_numbers, &env);

	if (flags & (PTM_CHECK_SC | PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO | PTM_CHECK_BCC)) {

		ptm::normalize_vertices(num_input_points, env.points, ch_points);
		ch.ok = false;

		if (flags & PTM_CHECK_SC)
			ret = match_general(&ptm::structure_sc, ch_points, env.points, &ch, &res);

		if (flags & (PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO))
			ret = match_fcc_hcp_ico(ch_points, env.points, flags, &ch, &res);

		if (flags & PTM_CHECK_BCC)
			ret = match_general(&ptm::structure_bcc, ch_points, env.points, &ch, &res);
	}

	if (flags & (PTM_CHECK_DCUB | PTM_CHECK_DHEX)) {

		const int num_inner = 4, num_outer = 3;

		ret = ptm::calculate_two_shell_neighbour_ordering(
				local_handle, num_input_points, num_inner, num_outer, &env, &dmn_env);
		if (ret == 0) {
			ptm::normalize_vertices(PTM_NUM_NBRS_DCUB + 1, dmn_env.points, ch_points);
			ch.ok = false;

			ret = match_dcub_dhex(ch_points, dmn_env.points, flags, &ch, &res);
		}
	}

	if (flags & PTM_CHECK_GRAPHENE) {

		const int num_inner = 3, num_outer = 2;

		ret = ptm::calculate_two_shell_neighbour_ordering(
				local_handle, num_input_points, num_inner, num_outer, &env, &grp_env);
		if (ret == 0) {
			ret = match_graphene(grp_env.points, &res);
		}
	}

	if (res.ref_struct == NULL)
		return PTM_NO_ERROR;

	ptm::atomicenv_t* res_env = &env;
	if (res.ref_struct->type == PTM_MATCH_DCUB || res.ref_struct->type == PTM_MATCH_DHEX)
		res_env = &dmn_env;
	else if (res.ref_struct->type == PTM_MATCH_GRAPHENE)
		res_env = &grp_env;

	output_data(	&res, res_env, output_conventional_orientation, p_type, p_alloy_type, p_scale,
			p_rmsd, q, F, F_res, U, P, p_interatomic_distance,
			p_lattice_constant, p_best_template_index, p_best_template, output_indices);

	return PTM_NO_ERROR;
}

