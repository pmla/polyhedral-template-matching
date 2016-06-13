#ifndef INDEX_PTM_HPP
#define INDEX_PTM_HPP

#include <stdint.h>
#include <stdbool.h>

//--------------------------------
//    definitions
//--------------------------------
#define PTM_NO_ERROR	0


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


//--------------------------------
//    internal definitions
//--------------------------------
#define MAX_NBRS	14
#define MAX_POINTS	15
#define MAX_FACETS	24

//--------------------------------
//    function declarations
//--------------------------------
#ifdef __cplusplus
extern "C" {
#endif

typedef struct ptm_local_handle* ptm_local_handle_t;
ptm_local_handle_t ptm_initialize_local();
void ptm_uninitialize_local(ptm_local_handle_t ptr);

int ptm_initialize_global();
int ptm_index(	ptm_local_handle_t local_handle, int num_points, double* atomic_positions, int32_t* atomic_numbers, int32_t flags, bool topological_ordering,						//inputs
		int32_t* p_type, int32_t* p_alloy_type, double* p_scale, double* p_rmsd, double* q, double* F, double* F_res, double* U, double* P, int8_t* mapping, double* p_lattice_constant);	//outputs

//--------------------------------
//    template structures
//--------------------------------

//these point sets have barycentre {0, 0, 0} and are scaled such that the mean neighbour distance is 1
const double ptm_template_sc[7][3] = {		{  0.            ,  0.            ,  0.             },
						{  0.            ,  0.            , -1.             },
						{  0.            ,  0.            ,  1.             },
						{  0.            , -1.            ,  0.             },
						{  0.            ,  1.            ,  0.             },
						{ -1.            ,  0.            ,  0.             },
						{  1.            ,  0.            ,  0.             }	};

const double ptm_template_fcc[13][3] = {	{  0.            ,  0.            ,  0.            },
						{  0.            ,  0.707106781187,  0.707106781187 },
						{  0.            , -0.707106781187, -0.707106781187 },
						{  0.            ,  0.707106781187, -0.707106781187 },
						{  0.            , -0.707106781187,  0.707106781187 },
						{  0.707106781187,  0.            ,  0.707106781187 },
						{ -0.707106781187,  0.            , -0.707106781187 },
						{  0.707106781187,  0.            , -0.707106781187 },
						{ -0.707106781187,  0.            ,  0.707106781187 },
						{  0.707106781187,  0.707106781187,  0.             },
						{ -0.707106781187, -0.707106781187,  0.             },
						{  0.707106781187, -0.707106781187,  0.             },
						{ -0.707106781187,  0.707106781187,  0.             }	};

const double ptm_template_hcp[13][3] = {	{  0.            ,  0.            ,  0.            },
						{  0.707106781186,  0.            ,  0.707106781186 },
						{ -0.235702260395, -0.942809041583, -0.235702260395 },
						{  0.707106781186,  0.707106781186,  0.             },
						{ -0.235702260395, -0.235702260395, -0.942809041583 },
						{  0.            ,  0.707106781186,  0.707106781186 },
						{ -0.942809041583, -0.235702260395, -0.235702260395 },
						{ -0.707106781186,  0.707106781186,  0.             },
						{  0.            ,  0.707106781186, -0.707106781186 },
						{  0.707106781186,  0.            , -0.707106781186 },
						{  0.707106781186, -0.707106781186,  0.             },
						{ -0.707106781186,  0.            ,  0.707106781186 },
						{  0.            , -0.707106781186,  0.707106781186 }	};

const double ptm_template_ico[13][3] = {	{  0.            ,  0.            ,  0.            },
						{  0.            ,  0.525731112119,  0.850650808352 },
						{  0.            , -0.525731112119, -0.850650808352 },
						{  0.            ,  0.525731112119, -0.850650808352 },
						{  0.            , -0.525731112119,  0.850650808352 },
						{ -0.525731112119, -0.850650808352,  0.             },
						{  0.525731112119,  0.850650808352,  0.             },
						{  0.525731112119, -0.850650808352,  0.             },
						{ -0.525731112119,  0.850650808352,  0.             },
						{ -0.850650808352,  0.            , -0.525731112119 },
						{  0.850650808352,  0.            ,  0.525731112119 },
						{  0.850650808352,  0.            , -0.525731112119 },
						{ -0.850650808352,  0.            ,  0.525731112119 }	};

const double ptm_template_bcc[15][3] = {	{  0.            ,  0.            ,  0.            },
						{ -0.541451884327, -0.541451884327, -0.541451884327 },
						{  0.541451884327,  0.541451884327,  0.541451884327 },
						{  0.541451884327, -0.541451884327, -0.541451884327 },
						{ -0.541451884327,  0.541451884327,  0.541451884327 },
						{ -0.541451884327,  0.541451884327, -0.541451884327 },
						{  0.541451884327, -0.541451884327,  0.541451884327 },
						{ -0.541451884327, -0.541451884327,  0.541451884327 },
						{  0.541451884327,  0.541451884327, -0.541451884327 },
						{  0.            ,  0.            , -1.082903768655 },
						{  0.            ,  0.            ,  1.082903768655 },
						{  0.            , -1.082903768655,  0.             },
						{  0.            ,  1.082903768655,  0.             },
						{ -1.082903768655,  0.            ,  0.             },
						{  1.082903768655,  0.            ,  0.             }	};

#ifdef __cplusplus
}
#endif

#endif

