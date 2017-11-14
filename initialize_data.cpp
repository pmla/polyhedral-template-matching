#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <algorithm>
#include "initialize_data.hpp"
#include "ptm_functions.h"


//refdata_t structure_sc =  { .type = PTM_MATCH_SC,  .num_nbrs =  6, .num_facets =  8, .max_degree = 4, .num_graphs = NUM_SC_GRAPHS,  .graphs = graphs_sc,  .points = ptm_template_sc,  .penrose = penrose_sc , .mapping = mapping_sc };
refdata_t structure_sc =  { PTM_MATCH_SC,   6,  8, 4, NUM_SC_GRAPHS,  graphs_sc,  ptm_template_sc,  penrose_sc , mapping_sc };
refdata_t structure_fcc = { PTM_MATCH_FCC, 12, 20, 6, NUM_FCC_GRAPHS, graphs_fcc, ptm_template_fcc, penrose_fcc, mapping_fcc};
refdata_t structure_hcp = { PTM_MATCH_HCP, 12, 20, 6, NUM_HCP_GRAPHS, graphs_hcp, ptm_template_hcp, penrose_hcp, mapping_hcp};
refdata_t structure_ico = { PTM_MATCH_ICO, 12, 20, 6, NUM_ICO_GRAPHS, graphs_ico, ptm_template_ico, penrose_ico, mapping_ico};
refdata_t structure_bcc = { PTM_MATCH_BCC, 14, 24, 8, NUM_BCC_GRAPHS, graphs_bcc, ptm_template_bcc, penrose_bcc, mapping_bcc};


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

ptm_local_handle_t ptm_initialize_local()
{
	assert(ptm_initialized);
	return (ptm_local_handle_t)voronoi_initialize_local();
}

void ptm_uninitialize_local(ptm_local_handle_t ptr)
{
	voronoi_uninitialize_local(ptr);
}

