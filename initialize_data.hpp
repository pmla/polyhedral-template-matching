#ifndef INITIALIZE_DATA_HPP
#define INITIALIZE_DATA_HPP


#include <cstdint>
#include "graph_data.hpp"
#include "deformation_gradient.hpp"
#include "fundamental_mappings.hpp"
#include "neighbour_ordering.hpp"
#include "canonical.hpp"
#include "convex_hull_incremental.hpp"
#include "ptm_constants.h"


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
} refdata_t;


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

