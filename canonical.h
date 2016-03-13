#ifndef CANONICAL_H
#define CANONICAL_H

#include <stdint.h>

uint64_t canonical_form(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree, int8_t* canonical_labelling);

#endif

