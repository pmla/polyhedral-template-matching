#ifndef PTM_CONSTANTS_H
#define PTM_CONSTANTS_H

//------------------------------------
//    definitions
//------------------------------------
#define PTM_NO_ERROR	0


#define PTM_CHECK_FCC	(1 << 0)
#define PTM_CHECK_HCP	(1 << 1)
#define PTM_CHECK_BCC	(1 << 2)
#define PTM_CHECK_ICO	(1 << 3)
#define PTM_CHECK_SC	(1 << 4)
#define PTM_CHECK_ALL	(PTM_CHECK_SC | PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO | PTM_CHECK_BCC)

#define PTM_MATCH_NONE	0
#define PTM_MATCH_FCC	1
#define PTM_MATCH_HCP	2
#define PTM_MATCH_BCC	3
#define PTM_MATCH_ICO	4
#define PTM_MATCH_SC	5

#define PTM_ALLOY_NONE		0
#define PTM_ALLOY_PURE		1
#define PTM_ALLOY_L10		2
#define PTM_ALLOY_L12_CU	3
#define PTM_ALLOY_L12_AU	4
#define PTM_ALLOY_B2		5


#define PTM_MAX_INPUT_POINTS 19
#define PTM_MAX_NBRS	14
#define PTM_MAX_POINTS	15
#define PTM_MAX_FACETS	24


//------------------------------------
//    number of neighbours
//------------------------------------
#define PTM_NUM_NBRS_FCC 12
#define PTM_NUM_NBRS_HCP 12
#define PTM_NUM_NBRS_BCC 14
#define PTM_NUM_NBRS_ICO 12
#define PTM_NUM_NBRS_SC  6

#define PTM_NUM_POINTS_FCC (PTM_NUM_NBRS_FCC + 1)
#define PTM_NUM_POINTS_HCP (PTM_NUM_NBRS_HCP + 1)
#define PTM_NUM_POINTS_BCC (PTM_NUM_NBRS_BCC + 1)
#define PTM_NUM_POINTS_ICO (PTM_NUM_NBRS_ICO + 1)
#define PTM_NUM_POINTS_SC  (PTM_NUM_NBRS_SC  + 1)

const int ptm_num_nbrs[6] = {0, PTM_NUM_NBRS_FCC, PTM_NUM_NBRS_HCP, PTM_NUM_NBRS_BCC, PTM_NUM_NBRS_ICO, PTM_NUM_NBRS_SC};

//------------------------------------
//    template structures
//------------------------------------

//these point sets have barycentre {0, 0, 0} and are scaled such that the mean neighbour distance is 1

const double ptm_template_fcc[PTM_NUM_POINTS_FCC][3] = {	{  0.            ,  0.            ,  0.             },
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

const double ptm_template_hcp[PTM_NUM_POINTS_HCP][3] = {	{  0.            ,  0.            ,  0.             },
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

const double ptm_template_bcc[PTM_NUM_POINTS_BCC][3] = {	{  0.            ,  0.            ,  0.             },
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

const double ptm_template_ico[PTM_NUM_POINTS_ICO][3] = {	{  0.            ,  0.            ,  0.             },
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

const double ptm_template_sc[PTM_NUM_POINTS_SC][3] = {		{  0.            ,  0.            ,  0.             },
								{  0.            ,  0.            , -1.             },
								{  0.            ,  0.            ,  1.             },
								{  0.            , -1.            ,  0.             },
								{  0.            ,  1.            ,  0.             },
								{ -1.            ,  0.            ,  0.             },
								{  1.            ,  0.            ,  0.             }	};
#endif

