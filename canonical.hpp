#ifndef CANONICAL_HPP
#define CANONICAL_HPP

#include <cstdint>

uint64_t canonical_form(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree, int8_t* canonical_labelling);

#endif

