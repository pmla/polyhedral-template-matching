#ifndef GRAPH_TOOLS_HPP
#define GRAPH_TOOLS_HPP

#include <cstdint>
#include "ptm_constants.h"

bool build_facet_map(int num_facets, int8_t facets[][3], int8_t common[PTM_MAX_NBRS][PTM_MAX_NBRS]);
int graph_degree(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree);

#endif

