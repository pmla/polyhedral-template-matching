#ifndef PTM_FUNCTIONS_H
#define PTM_FUNCTIONS_H

#include <stdint.h>
#include <stdbool.h>
#include "ptm_initialize_data.h"
#include "ptm_constants.h"


//------------------------------------
//    function declarations
//------------------------------------
#ifdef __cplusplus
extern "C" {
#endif


int ptm_index(	ptm_local_handle_t local_handle,
		size_t atom_index, int (get_neighbours)(void* vdata, size_t atom_index, int num, size_t* nbr_indices, int32_t* numbers, double (*nbr_pos)[3]), void* nbrlist,
		int32_t flags, bool output_conventional_orientation, //inputs
		int32_t* p_type, int32_t* p_alloy_type, double* p_scale, double* p_rmsd, double* q, double* F, double* F_res, double* U, double* P, double* p_interatomic_distance, double* p_lattice_constant,
		size_t* output_indices);	//outputs


#ifdef __cplusplus
}
#endif

#endif

