#ifndef PTM_FUNCTIONS_HPP
#define PTM_FUNCTIONS_HPP

#include <stdint.h>
#include <stdbool.h>
#include "initialize_data.hpp"
#include "ptm_constants.h"


//------------------------------------
//    function declarations
//------------------------------------
#ifdef __cplusplus
extern "C" {
#endif


int ptm_index(	ptm_local_handle_t local_handle, int num_points, double* atomic_positions, int32_t* atomic_numbers, int32_t flags, bool topological_ordering,										//inputs
		int32_t* p_type, int32_t* p_alloy_type, double* p_scale, double* p_rmsd, double* q, double* F, double* F_res, double* U, double* P, int8_t* mapping, double* p_interatomic_distance, double* p_lattice_constant);	//outputs


#ifdef __cplusplus
}
#endif

#endif

