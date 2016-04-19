#ifndef INDEX_PTM_H
#define INDEX_PTM_H

#include <stdint.h>
#include <stdbool.h>

#define PTM_CHECK_SC	(1 << 0)
#define PTM_CHECK_FCC	(1 << 1)
#define PTM_CHECK_HCP	(1 << 2)
#define PTM_CHECK_ICO	(1 << 3)
#define PTM_CHECK_BCC	(1 << 4)
#define PTM_CHECK_ALL	(PTM_CHECK_SC | PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO | PTM_CHECK_BCC)

#define PTM_MATCH_NONE	0
#define PTM_MATCH_SC	1
#define PTM_MATCH_FCC	2
#define PTM_MATCH_HCP	3
#define PTM_MATCH_ICO	4
#define PTM_MATCH_BCC	5

#define PTM_ALLOY_NONE		0
#define PTM_ALLOY_PURE		1
#define PTM_ALLOY_L10		2
#define PTM_ALLOY_L12_CU	3
#define PTM_ALLOY_L12_AU	4
#define PTM_ALLOY_B2		5


#ifdef __cplusplus
extern "C" {
#endif

int initialize_PTM();
void index_PTM(	int num_points, double* points, int32_t* numbers, int32_t flags, bool topological_ordering,						//inputs
		int32_t* p_type, int32_t* p_alloy_type, double* p_scale, double* p_rmsd, double* q, double* F, double* F_res, double* U, double* P);	//outputs

#ifdef __cplusplus
}
#endif

#endif

