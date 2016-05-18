#ifndef NEIGHBOUR_ORDERING_H
#define NEIGHBOUR_ORDERING_H

#include "index_ptm.h"

#ifdef __cplusplus
extern "C" {
#endif

int calculate_neighbour_ordering(ptm_local_handle_t, int num_points, const double (*_points)[3], int8_t* ordering);

#ifdef __cplusplus
}
#endif

#endif

