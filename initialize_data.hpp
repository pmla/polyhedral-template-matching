#ifndef INITIALIZE_DATA_HPP
#define INITIALIZE_DATA_HPP


#include <cstdint>
#include "graph_data.hpp"
#include "graph_tools.hpp"
#include "deformation_gradient.hpp"
#include "fundamental_mappings.hpp"
#include "neighbour_ordering.hpp"
#include "canonical_coloured.hpp"
#include "convex_hull_incremental.hpp"

#include <map>
#include <array>
#include <vector>
using namespace std;

typedef map< array<int8_t, 2 * PTM_MAX_EDGES>, std::vector<graph_t*> > graphmap_t;

extern graphmap_t map_sc;
extern graphmap_t map_fcc;
extern graphmap_t map_hcp;
extern graphmap_t map_ico;
extern graphmap_t map_bcc;
extern graphmap_t map_dcub;
extern graphmap_t map_dhex;

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
	const int8_t (*mapping)[PTM_MAX_POINTS];

	graphmap_t* graphmap;
} refdata_t;


//refdata_t structure_sc =  { .type = PTM_MATCH_SC,  .num_nbrs =  6, .num_facets =  8, .max_degree = 4, .num_graphs = NUM_SC_GRAPHS,  .graphs = graphs_sc,  .points = ptm_template_sc,  .penrose = penrose_sc , .mapping = mapping_sc };
const refdata_t structure_sc =   { PTM_MATCH_SC,    6,  8, 4, NUM_SC_GRAPHS,   graphs_sc,   ptm_template_sc,   penrose_sc,   mapping_sc,   &map_sc  };
const refdata_t structure_fcc =  { PTM_MATCH_FCC,  12, 20, 6, NUM_FCC_GRAPHS,  graphs_fcc,  ptm_template_fcc,  penrose_fcc,  mapping_fcc,  &map_fcc };
const refdata_t structure_hcp =  { PTM_MATCH_HCP,  12, 20, 6, NUM_HCP_GRAPHS,  graphs_hcp,  ptm_template_hcp,  penrose_hcp,  mapping_hcp,  &map_hcp };
const refdata_t structure_ico =  { PTM_MATCH_ICO,  12, 20, 6, NUM_ICO_GRAPHS,  graphs_ico,  ptm_template_ico,  penrose_ico,  mapping_ico,  &map_ico };
const refdata_t structure_bcc =  { PTM_MATCH_BCC,  14, 24, 8, NUM_BCC_GRAPHS,  graphs_bcc,  ptm_template_bcc,  penrose_bcc,  mapping_bcc,  &map_bcc };
const refdata_t structure_dcub = { PTM_MATCH_DCUB, 16, 28, 7, NUM_DCUB_GRAPHS, graphs_dcub, ptm_template_dcub, penrose_dcub, NULL, &map_dcub};
const refdata_t structure_dhex = { PTM_MATCH_DHEX, 16, 28, 7, NUM_DHEX_GRAPHS, graphs_dhex, ptm_template_dhex, penrose_dhex, NULL, &map_dhex};

const refdata_t structure_diac = { PTM_MATCH_DCUB, 12, 20, 6, NUM_FCC_GRAPHS, graphs_fcc, ptm_template_fcc, penrose_fcc, mapping_fcc, &map_dcub};
const refdata_t structure_diah = { PTM_MATCH_DHEX, 12, 20, 6, NUM_HCP_GRAPHS, graphs_hcp, ptm_template_hcp, penrose_hcp, mapping_hcp, &map_dhex};


#ifdef __cplusplus
extern "C" {
#endif

typedef struct ptm_local_handle* ptm_local_handle_t;
ptm_local_handle_t ptm_initialize_local();
void ptm_uninitialize_local(ptm_local_handle_t ptr);

int ptm_initialize_global();

//------------------------------------
//    global initialization switch
//------------------------------------
extern bool ptm_initialized;


#ifdef __cplusplus
}
#endif


#endif

