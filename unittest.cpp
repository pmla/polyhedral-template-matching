/*
	structures:
		+ no match for each of them
		nonzero rmsd

	add permutations of points
		with alloy numbers too

done:
	SC, FCC, HCP, ICO, BCC
		translation offset
		scale
		rotation
		zero rmsd
		orientations + fundamental zones
		strain rotation must equal to rmsd rotation if no anisotropic distortion, when mapped into fundamental zone

	alloys
		disordered + ordered

	deformation gradients

	strains + polar decomposition rotations

	check det(U) > 0
*/

#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "normalize_vertices.hpp"
#include "quat.hpp"
#include "ptm_functions.h"


#define RADIANS(x) (2.0 * M_PI * (x) / 360.0)
#define DEGREES(x) (360 * (x) / (2.0 * M_PI))


typedef struct
{
	int type;
	int check;
	int num_points;
	const double (*points)[3];
} structdata_t;

//                              { .type = PTM_MATCH_SC,  .check = PTM_CHECK_SC,  .num_points =  7, .points = ptm_template_sc },
structdata_t structdata[7] =  {	{ PTM_MATCH_FCC,  PTM_CHECK_FCC,  13, ptm_template_fcc  },
				{ PTM_MATCH_HCP,  PTM_CHECK_HCP,  13, ptm_template_hcp  },
				{ PTM_MATCH_BCC,  PTM_CHECK_BCC,  15, ptm_template_bcc  },
				{ PTM_MATCH_ICO,  PTM_CHECK_ICO,  13, ptm_template_ico  },
				{ PTM_MATCH_SC,   PTM_CHECK_SC,    7, ptm_template_sc   },
				{ PTM_MATCH_DCUB, PTM_CHECK_DCUB, 17, ptm_template_dcub },
				{ PTM_MATCH_DHEX, PTM_CHECK_DHEX, 17, ptm_template_dhex }, };

typedef struct
{
	int32_t type;
	int32_t numbers[PTM_MAX_POINTS];
} alloytest_t;

alloytest_t sc_alloy_tests[] = {

	{ PTM_ALLOY_NONE,   {-1}},	//no test
};

alloytest_t fcc_alloy_tests[] = {

	{ PTM_ALLOY_NONE,   {-1}},	//no test

	{ PTM_ALLOY_NONE,   {0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},	//pure -defect
	{ PTM_ALLOY_NONE,   {4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 4, 4, 4}},	//pure -defect
	{ PTM_ALLOY_NONE,   {3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 7}},	//L10 -defect
	{ PTM_ALLOY_NONE,   {3, 0, 3, 0, 0, 3, 3, 3, 3, 0, 0, 0, 0}},	//L10 -defect
	{ PTM_ALLOY_NONE,   {3, 0, 1, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3}},	//L10 -defect
	{ PTM_ALLOY_NONE,   {0, 3, 1, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0}},	//L12_CU -defect
	{ PTM_ALLOY_NONE,   {5, 5, 3, 5, 5, 3, 3, 3, 3, 5, 5, 5, 5}},	//L12_CU -defect
	{ PTM_ALLOY_NONE,   {9, 9, 2, 9, 9, 9, 9, 9, 9, 3, 3, 3, 3}},	//L12_CU -defect
	{ PTM_ALLOY_NONE,   {4, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},	//L12_AU -defect
	{ PTM_ALLOY_NONE,   {1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}},	//L12_AU -defect

	{ PTM_ALLOY_PURE,   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
	{ PTM_ALLOY_PURE,   {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},
	{ PTM_ALLOY_L10,    {3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0}},
	{ PTM_ALLOY_L10,    {3, 0, 0, 0, 0, 3, 3, 3, 3, 0, 0, 0, 0}},
	{ PTM_ALLOY_L10,    {3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3}},
	{ PTM_ALLOY_L12_CU, {0, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0}},
	{ PTM_ALLOY_L12_CU, {5, 5, 5, 5, 5, 3, 3, 3, 3, 5, 5, 5, 5}},
	{ PTM_ALLOY_L12_CU, {9, 9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 3, 3}},
	{ PTM_ALLOY_L12_AU, {4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
	{ PTM_ALLOY_L12_AU, {1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}},
};

alloytest_t hcp_alloy_tests[] = {

	{ PTM_ALLOY_NONE,   {-1}},	//no test
};

alloytest_t ico_alloy_tests[] = {

	{ PTM_ALLOY_NONE,   {-1}},	//no test
};

alloytest_t bcc_alloy_tests[] = {

	{ PTM_ALLOY_NONE,   {-1}},	//no test

	{ PTM_ALLOY_NONE,   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}},	//pure -defect
	{ PTM_ALLOY_NONE,   {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},	//pure -defect
	{ PTM_ALLOY_NONE,   {4, 0, 0, 0, 0, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4}},	//B2 -defect
	{ PTM_ALLOY_NONE,   {1, 2, 2, 2, 2, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1}},	//B2 -defect

	{ PTM_ALLOY_PURE,   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
	{ PTM_ALLOY_PURE,   {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},
	{ PTM_ALLOY_B2,     {4, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4}},
	{ PTM_ALLOY_B2,     {1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1}},
};

alloytest_t dcub_alloy_tests[] = {

	{ PTM_ALLOY_NONE,   {-1}},	//no test

	{ PTM_ALLOY_NONE,   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}},	//pure -defect
	{ PTM_ALLOY_NONE,   {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},	//pure -defect
	{ PTM_ALLOY_NONE,   {4, 0, 0, 0, 0, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4}},	//SiC -defect
	{ PTM_ALLOY_NONE,   {1, 2, 2, 2, 2, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1}},	//SiC -defect

	{ PTM_ALLOY_PURE,   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
	{ PTM_ALLOY_PURE,   {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},
	{ PTM_ALLOY_SIC,    {4, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},
	{ PTM_ALLOY_SIC,    {1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},

};


alloytest_t dhex_alloy_tests[] = {

	{ PTM_ALLOY_NONE,   {-1}},	//no test
};


typedef struct
{
	bool strain;
	double pre[9];
	double post[9];
} quattest_t;

quattest_t cubic_qtest[] = {	{	false, {1.00, 0.00, -0.00, 0.00},			{1.00, 0.00, -0.00, 0.00}			},
				{	false, {0.00, 1.00, -0.00, 0.00},			{1.00, 0.00, -0.00, 0.00}			},
				{	false, {0.00, 0.00, -1.00, 0.00},			{1.00, 0.00, -0.00, 0.00}			},
				{	false, {0.00, 0.00, -0.00, 1.00},			{1.00, 0.00, -0.00, 0.00}			},
				{	false, {0.987070, 0.020780, -0.031171, 0.155853},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.683270, 0.712658, 0.088164, 0.132246},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.720005, -0.095511, 0.675923, 0.124899},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.587759, -0.007347, -0.036735, 0.808168},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.808168, 0.036735, -0.007347, -0.587759},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.675923, 0.124899, -0.720005, 0.095511},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.712658, -0.683270, -0.132246, 0.088164},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.420803, 0.410413, 0.545486, 0.597437},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.576656, 0.441584, 0.566266, -0.389633},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.389633, 0.566266, -0.441584, 0.576656},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.545486, 0.597437, -0.420803, -0.410413},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.441584, -0.576656, 0.389633, 0.566266},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.597437, -0.545486, 0.410413, -0.420803},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.410413, -0.420803, -0.597437, 0.545486},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.566266, -0.389633, -0.576656, -0.441584},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.020780, -0.987070, -0.155853, -0.031171},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.007347, 0.587759, 0.808168, 0.036735},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.124899, -0.675923, -0.095511, -0.720005},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.095511, 0.720005, 0.124899, -0.675923},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.036735, -0.808168, 0.587759, -0.007347},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.031171, -0.155853, 0.987070, 0.020780},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.088164, 0.132246, -0.683270, -0.712658},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.132246, -0.088164, 0.712658, -0.683270},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	false, {0.155853, 0.031171, 0.020780, -0.987070},	{0.987070, 0.020780, -0.031171, 0.155853}	},
				{	true,  {1.1, 0.0, 0.0, -0.0, 1.0, -0.0, -0.0, 0.0, 1.0},	{-1, -1, -1, -1}			},
				{	true,  {1.05, 0.0, 0.0, -0.05, 0.95, -0.0, -0.02, 0.06, 0.94},	{-1, -1, -1, -1}			}	};

quattest_t ico_qtest[] = {	{	false, {0.960379, 0.109103, 0.249483, -0.059392},	{0.960379, 0.109103, 0.249483, -0.059392}	},
				{	false, {0.711212, 0.456880, 0.518044, -0.130648},	{0.975588, 0.084523, 0.103873, -0.174053}	},
				{	false, {0.936649, 0.067820, 0.187887, 0.287729},	{0.959690, 0.049837, -0.103527, -0.256504}	},
				{	false, {0.638434, 0.650457, 0.048934, 0.408550},	{0.967981, 0.191893, 0.036316, 0.157704}	},
				{	false, {0.821165, -0.149986, 0.542985, -0.091435},	{0.982177, 0.086695, 0.000446, -0.166771}	},
				{	false, {0.639462, -0.154114, -0.578165, 0.482766},	{0.937381, 0.015219, -0.347198, 0.023211}	},
				{	false, {0.753493, 0.595008, 0.139240, -0.242540},	{0.982042, 0.147653, 0.050050, 0.106243}	},
				{	false, {0.701011, -0.659762, 0.158265, 0.219660},	{0.964889, -0.232160, 0.033991, -0.118048}	},
				{	false, {0.769565, 0.183283, 0.380132, -0.479246},	{0.979681, 0.190250, -0.021917, -0.059573}	},
				{	false, {0.931242, 0.075646, -0.356282, 0.011355},	{0.954908, -0.232248, 0.173874, -0.063088}	},
				{	false, {0.867829, 0.097112, -0.422785, -0.242270},	{0.953871, -0.057962, -0.122423, 0.267923}	},
				{	false, {0.920769, -0.141468, -0.337810, -0.134371},	{0.957539, 0.237269, 0.145568, -0.075053}	},
				{	false, {0.697680, -0.476591, -0.029266, -0.534085},	{0.967772, -0.045774, -0.143445, -0.201856}	},
				{	false, {0.407309, -0.819150, 0.372159, -0.156812},	{0.981364, -0.128513, -0.066650, -0.126357}	},
				{	false, {0.503983, -0.834038, -0.158844, -0.158590},	{0.975828, 0.039719, -0.051984, -0.208520}	},
				{	false, {0.579265, -0.376209, -0.513427, -0.509227},	{0.989064, 0.099428, -0.033590, 0.103628}	},
				{	false, {0.612812, -0.641471, -0.339335, 0.312774},	{0.953196, -0.001049, -0.027610, -0.301087}	},
				{	false, {0.436573, -0.615735, 0.380584, -0.534257},	{0.983574, -0.166417, 0.012745, 0.068733}	},
				{	false, {0.269644, -0.465014, 0.604212, 0.588202},	{0.963536, -0.105690, 0.228878, 0.089680}	},
				{	false, {0.684766, -0.095996, -0.076056, -0.718399},	{0.953244, 0.102076, -0.182363, 0.218290}	},
				{	false, {0.451921, -0.281443, -0.099526, 0.840626},	{0.993012, 0.079448, -0.017687, 0.085457}	},
				{	false, {0.185681, -0.038838, -0.901638, 0.388668},	{0.942386, -0.055237, -0.312601, 0.105535}	},
				{	false, {0.442598, 0.081466, -0.678884, -0.580161},	{0.949807, 0.300306, -0.006548, -0.087403}	},
				{	false, {0.610532, -0.157028, 0.707014, 0.320507},	{0.976294, -0.037697, -0.188948, 0.098627}	},
				{	false, {0.517194, -0.052227, 0.820043, -0.239400},	{0.996004, 0.033614, 0.007742, 0.082374}	},
				{	false, {0.527468, 0.290749, -0.094791, -0.792627},	{0.994829, -0.094310, -0.037682, 0.001125}	},
				{	false, {0.518144, 0.525232, 0.182485, 0.649891},	{0.947150, -0.045132, 0.315337, -0.037851}	},
				{	false, {0.425858, 0.546800, -0.225130, -0.684815},	{0.941301, 0.290314, 0.169372, 0.031357}	},
				{	false, {0.729071, 0.460371, -0.387695, 0.325893},	{0.951515, -0.103449, 0.237927, -0.165251}	},
				{	false, {0.585124, 0.389119, 0.542312, -0.460558},	{0.988557, -0.057126, 0.014314, 0.138880}	},
				{	false, {0.552261, 0.398459, 0.402666, 0.611636},	{0.982511, 0.027584, -0.181387, 0.031791}	},
				{	false, {0.234592, 0.919558, 0.314282, 0.024627},	{0.958352, 0.277601, 0.064725, -0.017586}	},
				{	false, {0.686597, 0.680007, -0.257241, -0.001496},	{0.972928, -0.215003, 0.084760, 0.001273}	},
				{	false, {0.301914, 0.927941, 0.106506, 0.190868},	{0.939450, -0.010757, 0.342467, -0.005811}	},
				{	false, {0.184193, -0.202322, 0.518320, 0.810237},	{0.971574, -0.076732, -0.095608, 0.202523}	},
				{	false, {0.454091, 0.354246, 0.784990, -0.228259},	{0.952515, -0.302243, -0.010663, 0.035368}	},
				{	false, {0.176659, 0.156611, -0.535053, 0.811162},	{0.978361, 0.075681, 0.049691, 0.186048}	},
				{	false, {0.389189, 0.795539, 0.281348, -0.369448},	{0.948595, 0.071649, -0.011938, 0.308044}	},
				{	false, {0.022236, -0.390218, 0.918085, -0.065986},	{0.944727, -0.162850, 0.232722, -0.163740}	},
				{	false, {0.389599, -0.780594, -0.283174, 0.398370},	{0.951091, 0.215562, -0.155514, 0.157396}	},
				{	false, {0.518603, 0.417864, -0.716708, -0.206810},	{0.949018, 0.037138, 0.301489, -0.084202}	},
				{	false, {0.021515, -0.148317, 0.553731, -0.819098},	{0.946178, -0.007404, 0.280346, -0.161551}	},
				{	false, {0.421708, -0.763995, 0.284758, -0.396721},	{0.946760, 0.247461, 0.149039, -0.142114}	},
				{	false, {0.158295, -0.091623, -0.567854, -0.802552},	{0.982121, -0.086440, -0.022205, -0.165750}	},
				{	false, {0.074503, -0.323527, -0.932109, 0.144753},	{0.938878, -0.179832, -0.155387, 0.249047}	},
				{	false, {0.160104, 0.934696, 0.317340, 0.002102},	{0.934696, -0.160104, -0.002102, 0.317340}	},
				{	false, {0.008704, 0.835398, -0.501689, 0.224370},	{0.996030, 0.035803, 0.080985, 0.009134}	},
				{	false, {0.127408, -0.909278, 0.297157, 0.262064},	{0.965183, 0.142280, -0.132672, 0.174863}	},
				{	false, {0.173455, -0.696049, -0.643370, -0.267404},	{0.967433, 0.075218, 0.085484, 0.226073}	},
				{	false, {0.103343, 0.900503, 0.366757, -0.209532},	{0.976650, -0.075038, -0.160427, -0.121604}	},
				{	false, {0.135899, -0.426302, -0.059631, 0.892324},	{0.953483, -0.159551, 0.143271, -0.211864}	},
				{	false, {0.122320, 0.475403, 0.066330, 0.868695},	{0.960987, 0.153620, -0.087538, -0.212701}	},
				{	false, {0.020906, 0.617289, -0.517511, -0.592199},	{0.947663, -0.246129, -0.196838, -0.051089}	},
				{	false, {0.038143, -0.586443, 0.281833, -0.758420},	{0.993887, 0.012715, 0.083446, 0.071162}	},
				{	false, {0.161701, 0.368764, -0.670972, 0.622626},	{0.968095, -0.218198, 0.122798, 0.010144}	},
				{	false, {0.046935, 0.321773, -0.793677, -0.514136},	{0.998600, 0.004602, 0.035962, 0.038528}	},
				{	false, {0.071151, -0.476847, -0.779461, -0.399994},	{0.977948, -0.044142, 0.172381, -0.109336}	},
				{	false, {0.113851, -0.380862, -0.616422, 0.679710},	{0.956244, -0.206504, 0.111718, -0.174564}	},
				{	false, {0.001683, -0.124987, 0.988455, -0.085627},	{0.988455, -0.085627, -0.001683, 0.124987}	},
				{	false, {0.220005, -0.085991, 0.026420, -0.971342},	{0.971342, 0.026420, 0.085991, 0.220005}	},
				{	true,  {1.1, 0.0, 0.0, -0.0, 1.0, -0.0, -0.0, 0.0, 1.0}, 	{-1, -1, -1, -1}			},
				{	true,  {1.05, 0.0, 0.0, -0.05, 0.95, -0.0, -0.02, 0.06, 0.94},	{-1, -1, -1, -1}			}	};


quattest_t hcp_qtest[] = {	{	false, {1.00, 0.00, -0.00, 0.00},			{1.00, 0.00, -0.00, 0.00}			},
				{	false, {0.813035, -0.166389, 0.181571, -0.527562},	{0.813035, -0.166389, 0.181571, -0.527562}	},
				{	false, {0.023227, -0.154944, -0.861953, -0.482172},	{0.761148, -0.255749, -0.582977, 0.124032}	},
				{	false, {0.631878, 0.693281, 0.293658, -0.184001},	{0.717408, -0.208128, 0.269531, -0.607751}	},
				{	false, {0.167480, -0.812443, 0.552710, 0.079998},	{0.921658, -0.056237, -0.334734, -0.187981}	},
				{	false, {0.002364, -0.201656, 0.941752, 0.269132},	{0.741391, 0.605179, 0.190268, -0.218852}	},
				{	false, {0.122621, -0.141644, 0.064344, 0.980184},	{0.831875, -0.402635, 0.334567, -0.184214}	},
				{	true,  {1.1, 0.0, 0.0, -0.0, 1.0, -0.0, -0.0, 0.0, 1.0}, 	{-1, -1, -1, -1}			},
				{	true,  {1.05, 0.0, 0.0, -0.05, 0.95, -0.0, -0.02, 0.06, 0.94},	{-1, -1, -1, -1}			}	};

//#ifdef DEBUG
#define ERROR(msg, code) print_error(__FILE__, __PRETTY_FUNCTION__, __LINE__, msg, code)
#define CLEANUP(msg, code) {ret = code; print_error(__FILE__, __PRETTY_FUNCTION__, __LINE__, msg, code); goto cleanup;}
//#else
//#define ERROR(msg, code) code
//#define CLEANUP(msg, code) {ret = code; goto cleanup;}
//#endif

static int print_error(const char* file, const char* function, int line, const char* msg, int error_code)
{
	printf("\n\nerror\tfile: %s\n", file);
	printf("\tline: %d\n", line);
	printf("\tfunction: %s\n", function);
	printf("\terror message: %s\n", msg);
	printf("\terror code: %d\n", error_code);
	(void)fflush(stdout);
	return error_code;
}

static void matvec(double* A, double* x, double* b)
{
	b[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
	b[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
	b[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}

static void matmul(double* A, double* x, double* b)
{
	b[0] = A[0] * x[0] + A[1] * x[3] + A[2] * x[6];
	b[3] = A[3] * x[0] + A[4] * x[3] + A[5] * x[6];
	b[6] = A[6] * x[0] + A[7] * x[3] + A[8] * x[6];

	b[1] = A[0] * x[1] + A[1] * x[4] + A[2] * x[7];
	b[4] = A[3] * x[1] + A[4] * x[4] + A[5] * x[7];
	b[7] = A[6] * x[1] + A[7] * x[4] + A[8] * x[7];

	b[2] = A[0] * x[2] + A[1] * x[5] + A[2] * x[8];
	b[5] = A[3] * x[2] + A[4] * x[5] + A[5] * x[8];
	b[8] = A[6] * x[2] + A[7] * x[5] + A[8] * x[8];
}

static bool check_matrix_equality(double* A, double* B, double tolerance)
{
	for (int i=0;i<9;i++)
		if (fabs(A[i] - B[i]) > tolerance)
			return false;
	return true;
}

static double matrix_determinant(double* A)
{
	return    A[0] * (A[4] * A[8] - A[5] * A[7])
		- A[1] * (A[3] * A[8] - A[5] * A[6])
		+ A[2] * (A[3] * A[7] - A[4] * A[6]);
}

static double nearest_neighbour_rmsd(int num, double scale, double* A, double (*input_points)[3], const double (*template_points)[3])
{
	//transform template
	double transformed_template[PTM_MAX_POINTS][3];
	for (int i=0;i<num;i++)
	{
		double row[3] = {0, 0, 0};
		matvec(A, (double*)template_points[i], row);
		memcpy(transformed_template[i], row, 3 * sizeof(double));
	}

	//translate and scale input points
	double points[PTM_MAX_POINTS][3];
	subtract_barycentre(num, input_points, points);
	for (int i=0;i<num;i++)
		for (int j=0;j<3;j++)
			points[i][j] *= scale;

	double acc = 0;
	for (int i=0;i<num;i++)
	{
		double x0 = points[i][0];
		double y0 = points[i][1];
		double z0 = points[i][2];

		double min_dist = INFINITY;
		for (int j=0;j<num;j++)
		{
			double x1 = transformed_template[j][0];
			double y1 = transformed_template[j][1];
			double z1 = transformed_template[j][2];

			double dx = x1 - x0;
			double dy = y1 - y0;
			double dz = z1 - z0;
			double dist = dx*dx + dy*dy + dz*dz;
			min_dist = std::min(min_dist, dist);
		}

		acc += min_dist;
	}

	return sqrt(fabs(acc / num));
}

static double mapped_neighbour_rmsd(int num, double scale, double* A, double (*input_points)[3], const double (*template_points)[3], int8_t* mapping)
{
	//transform template
	double transformed_template[PTM_MAX_POINTS][3];
	for (int i=0;i<num;i++)
	{
		double row[3] = {0, 0, 0};
		matvec(A, (double*)template_points[i], row);
		memcpy(transformed_template[i], row, 3 * sizeof(double));
	}

	//translate and scale input points
	double points[PTM_MAX_POINTS][3];
	subtract_barycentre(num, input_points, points);
	for (int i=0;i<num;i++)
		for (int j=0;j<3;j++)
			points[i][j] *= scale;

	double acc = 0;
	for (int i=0;i<num;i++)
	{
		double x0 = points[mapping[i]][0];
		double y0 = points[mapping[i]][1];
		double z0 = points[mapping[i]][2];

		double x1 = transformed_template[i][0];
		double y1 = transformed_template[i][1];
		double z1 = transformed_template[i][2];

		double dx = x1 - x0;
		double dy = y1 - y0;
		double dz = z1 - z0;
		double dist = dx*dx + dy*dy + dz*dz;
		acc += dist;
	}

	return sqrt(fabs(acc / num));
}

typedef double points_t[3];

uint64_t run_tests()
{
	int ret = 0;
	const double tolerance = 1E-5;
	double identity_matrix[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	int num_structures = sizeof(structdata) / sizeof(structdata_t);

	int num_alloy_tests[] = {	sizeof(fcc_alloy_tests) / sizeof(alloytest_t),
					sizeof(hcp_alloy_tests) / sizeof(alloytest_t),
					sizeof(bcc_alloy_tests) / sizeof(alloytest_t),
					sizeof(ico_alloy_tests) / sizeof(alloytest_t),
					sizeof(sc_alloy_tests) / sizeof(alloytest_t),
					sizeof(dcub_alloy_tests) / sizeof(alloytest_t),
					sizeof(dhex_alloy_tests) / sizeof(alloytest_t)	};

	alloytest_t* alloy_test[] = {	fcc_alloy_tests,
					hcp_alloy_tests,
					bcc_alloy_tests,
					ico_alloy_tests,
					sc_alloy_tests,
					dcub_alloy_tests,
					dhex_alloy_tests  };

	int num_quat_tests[] = {	sizeof(cubic_qtest) / sizeof(quattest_t),
					sizeof(hcp_qtest) / sizeof(quattest_t),
					sizeof(cubic_qtest) / sizeof(quattest_t),
					sizeof(ico_qtest) / sizeof(quattest_t),
					sizeof(cubic_qtest) / sizeof(quattest_t),
					1,
					1	};

	quattest_t* quat_test[] = {	cubic_qtest,
					hcp_qtest,
					cubic_qtest,
					ico_qtest,
					cubic_qtest,
					cubic_qtest,
					hcp_qtest	};
	int num_tests = 0;
	ptm_local_handle_t local_handle = ptm_initialize_local();

	//rotation matrix => quaternion
	{
		double U[9] = {	 0,  0,  1,
				 0,  1,  0,
				-1,  0,  0  };
		double q[4];
		rotation_matrix_to_quaternion(U, q);
		double ans[4] = {1 / sqrt(2), 0, 1 / sqrt(2), 0};
		for (int i = 0;i<4;i++)
			if (fabs(q[i] - ans[i]) >= tolerance)
				CLEANUP("failed on rotation matrix => quaternion conversion", -1)
	}

	//quaternion => rotation matrix
	{
		double q[4] = {1 / sqrt(2), 1 / sqrt(2), 0, 0};
		double U[9];
		quaternion_to_rotation_matrix(q, U);

		double ans[9] = { 1,  0,  0,
				  0,  0, -1,
				  0,  1,  0	};

		for (int i = 0;i<9;i++)
			if (fabs(U[i] - ans[i]) >= tolerance)
				CLEANUP("failed on quaternion => rotation matrix conversion", -1)
	}

/*
	//for (int it = 0;it<num_structures;it++)
int it = 2;
	{
		structdata_t* s = &structdata[it];
		double points[PTM_MAX_POINTS][3];
		double rotated_points[PTM_MAX_POINTS][3];

		memcpy(points, s->points, 3 * sizeof(double) * s->num_points);

		//for (int k=0;k<24;k++)
		for (int k=0;k<6;k++)
		{
			double U[9];
			//extern double generator_cubic[24][4];
			//extern double generator_icosahedral[60][4];
			extern double generator_hcp[6][4];
			double qtemp[4];
			//memcpy(qtemp, generator_cubic[k], 4 * sizeof(double));
			//memcpy(qtemp, generator_icosahedral[k], 4 * sizeof(double));
			memcpy(qtemp, generator_hcp[k], 4 * sizeof(double));
			qtemp[1] = -qtemp[1];
			qtemp[2] = -qtemp[2];
			qtemp[3] = -qtemp[3];
			quaternion_to_rotation_matrix(qtemp, U);

			for (int i=0;i<s->num_points;i++)
			{
				double row[3] = {0, 0, 0};
				matvec(U, (double*)points[i], row);
				memcpy(rotated_points[i], row, 3 * sizeof(double));
			}

			for (int i=0;i<s->num_points;i++)
			{
				double x0 = rotated_points[i][0];
				double y0 = rotated_points[i][1];
				double z0 = rotated_points[i][2];

				int bi = -1;
				double min_dist = INFINITY;
				for (int j=0;j<s->num_points;j++)
				{
					double x1 = points[j][0];
					double y1 = points[j][1];
					double z1 = points[j][2];

					double dx = x1 - x0;
					double dy = y1 - y0;
					double dz = z1 - z0;
					double dist = dx*dx + dy*dy + dz*dz;
					if (dist < min_dist)
						bi = j;
					min_dist = MIN(min_dist, dist);
				}

				printf("%d ", bi);
			}
			printf("\n");
		}
	}
exit(3);*/


	for (int it = 0;it<num_structures;it++)
	{
		for (int iq=0;iq<num_quat_tests[it];iq++)
		{
			quattest_t* qtest = quat_test[it];

			double qpre[4], qpost[4], rot[9];
			memcpy(qpre, qtest[iq].pre, 4 * sizeof(double));
			memcpy(qpost, qtest[iq].post, 4 * sizeof(double));

			if (qtest[iq].strain)
			{
				memcpy(rot, qtest[iq].pre, 9 * sizeof(double));
			}
			else
			{
				normalize_quaternion(qpre);
				normalize_quaternion(qpost);
				quaternion_to_rotation_matrix(qpre, rot);
			}


			structdata_t* s = &structdata[it];
			double scaled_only[PTM_MAX_POINTS][3] = {0};
			double points[PTM_MAX_POINTS][3];
			memcpy(points, s->points, 3 * sizeof(double) * s->num_points);

			for (int i=0;i<s->num_points;i++)
			{
				double row[3] = {0, 0, 0};
				matvec(rot, points[i], row);
				memcpy(points[i], row, 3 * sizeof(double));
			}

			double _x = s->points[1][0];
			double _y = s->points[1][1];
			double _z = s->points[1][2];
			double norm = sqrt(_x*_x + _y*_y + _z*_z) / 2;
			double rescale = 1 / norm;
			double offset[3] = {12.0, 45.6, 789.10};
			for (int i = 0;i<s->num_points;i++)
				for (int j = 0;j<3;j++)
					scaled_only[i][j] = points[i][j] * rescale;

			for (int i = 0;i<s->num_points;i++)
				for (int j = 0;j<3;j++)
					points[i][j] = points[i][j] * rescale + offset[j];


			int tocheck = 0;
			for (int i = 0;i<it+1;i++)
				if (structdata[i].num_points <= structdata[it].num_points)
					tocheck |= structdata[i].check;

			for (int ia=0;ia<num_alloy_tests[it];ia++)
			{
				int32_t* numbers = alloy_test[it][ia].numbers;
				if (numbers[0] == -1)
					numbers = NULL;

				for (int itop=0;itop<=1;itop++)
				{
					int8_t mapping[PTM_MAX_POINTS];
					bool topological = itop == 1;
					int32_t type, alloy_type;
					double scale, rmsd, interatomic_distance, lattice_constant;
					double q[4], F[9], F_res[3], U[9], P[9];
					ret = ptm_index(local_handle, tocheck, s->num_points, points, numbers, topological, &type, &alloy_type, &scale, &rmsd, q, F, F_res, U, P, mapping, &interatomic_distance, &lattice_constant);
					if (ret != PTM_NO_ERROR)
						CLEANUP("indexing failed", ret);

					num_tests++;

#ifdef DEBUG
					printf("type:\t\t%d\t(should be: %d)\n", type, s->type);
					printf("alloy type:\t%d\n", alloy_type);
					printf("scale:\t\t%f\n", scale);
					printf("rmsd:\t\t%f\n", rmsd);
					printf("quat: \t\t%.4f %.4f %.4f %.4f\n", q[0], q[1], q[2], q[3]);
					printf("qpost:\t\t%.4f %.4f %.4f %.4f\n", qpost[0], qpost[1], qpost[2], qpost[3]);
					printf("rot: %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", rot[0], rot[1], rot[2], rot[3], rot[4], rot[5], rot[6], rot[7], rot[8]);
					printf("U:   %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7], U[8]);
					printf("P:   %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8]);
					printf("F:   %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", F[0], F[1], F[2], F[3], F[4], F[5], F[6], F[7], F[8]);
					//printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", rot[0] - U[0], rot[1] - U[1], rot[2] - U[2], rot[3] - U[3], rot[4] - U[4], rot[5] - U[5], rot[6] - U[6], rot[7] - U[7], rot[8] - U[8]);
					printf("interatomic distance:\t\t%f\n", interatomic_distance);
#endif

					//check type
					if (type != s->type)
						CLEANUP("failed on type", -1);

					//check alloy type
					if (alloy_type != alloy_test[it][ia].type)
						CLEANUP("failed on alloy type", -1);

					//check U-matrix is right handed
					if (matrix_determinant(U) <= 0)
						CLEANUP("failed on U-matrix right-handedness test", -1);

					//check strain tensor is symmetric
					if (fabs(P[1] - P[3]) > tolerance)	CLEANUP("failed on strain tensor symmetry test", -1);
					if (fabs(P[2] - P[6]) > tolerance)	CLEANUP("failed on strain tensor symmetry test", -1);
					if (fabs(P[5] - P[7]) > tolerance)	CLEANUP("failed on strain tensor symmetry test", -1);

					//check polar decomposition
					double _F[9];
					matmul(P, U, _F);
					if (!check_matrix_equality(_F, F, tolerance))
						CLEANUP("failed on polar decomposition check", -1);

					double A[9];
					if (!qtest[iq].strain)
					{
						//check rmsd
						if (rmsd > tolerance)
							CLEANUP("failed on rmsd", -1);

						//check scale
						if (fabs(scale - 1 / rescale) > tolerance)
							CLEANUP("failed on scale", -1);

						//check deformation gradient equal to polar decomposition rotation
						if (!check_matrix_equality(F, U, tolerance))
							CLEANUP("failed on deformation gradient check", -1);

						//check strain tensor is identity
						if (!check_matrix_equality(P, identity_matrix, tolerance))
							CLEANUP("failed on P identity matrix check", -1);

						//check rotation
						if (quat_misorientation(q, qpost) > tolerance)
							CLEANUP("failed on disorientation", -1);

						//check deformation gradient disorientation
						double qu[4];
						rotation_matrix_to_quaternion(U, qu);
						if (quat_misorientation(qu, qpost) > tolerance)
							CLEANUP("failed on deformation gradient disorientation", -1);

						quaternion_to_rotation_matrix(q, A);

						double x = scaled_only[1][0];
						double y = scaled_only[1][1];
						double z = scaled_only[1][2];
						double iad = sqrt(x*x + y*y + z*z);
						if (fabs(iad - interatomic_distance) > tolerance)
							CLEANUP("failed on interatomic distance", -1);
					}
					else
					{
						memcpy(A, F, 9 * sizeof(double));
					}

					//check nearest neighbour rmsd
					double rmsd_approx = nearest_neighbour_rmsd(s->num_points, scale, A, points, s->points);
					if (fabs(rmsd_approx) > tolerance)
						CLEANUP("failed on rmsd nearest neighbour", -1);

					//check mapped neighbour rmsd
					double rmsd_mapped = mapped_neighbour_rmsd(s->num_points, scale, A, points, s->points, mapping);
					if (fabs(rmsd_mapped) > tolerance)
						CLEANUP("failed on rmsd mapped neighbour", -1);
				}
			}
		}
	}


	{
		double lc_points_sc[7][3] = {{0,0,0},{2,0,0},{-2,0,0},{0,2,0},{0,-2,0},{0,0,2},{0,0,-2}};
		double lc_points_fcc[13][3] = {{0,0,0},{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1},{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0}};
		double lc_points_bcc[15][3] = {{0,0,0},{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1},{2,0,0},{-2,0,0},{0,2,0},{0,-2,0},{0,0,2},{0,0,-2}};
		double lc_points_dcub[17][3] = {{0,0,0},{-0.5,0.5,0.5},{-0.5,-0.5,-0.5},{0.5,-0.5,0.5},{0.5,0.5,-0.5},{-1,0,1},{-1,1,0},{0,1,1},{-1,-1,0},{-1,0,-1},{0,-1,-1},{0,-1,1},{1,-1,0},{1,0,1},{0,1,-1},{1,0,-1},{1,1,0}};
		points_t* pdata[4] = {lc_points_fcc, lc_points_bcc, lc_points_sc, lc_points_dcub};

		//double* (*pdata)[3] = {lc_points_fcc, lc_points_bcc, lc_points_sc};

		int lcdat[4] = {0, 2, 4, 5};
		for (int i=0;i<4;i++)
		{
			structdata_t* s = &structdata[lcdat[i]];

			int32_t type;
			double scale, rmsd, interatomic_distance, lattice_constant, q[4];
			ret = ptm_index(local_handle, s->check, s->num_points, pdata[i], NULL, false, &type, NULL, &scale, &rmsd, q, NULL, NULL,  NULL, NULL, NULL, &interatomic_distance, &lattice_constant);
			if (ret != PTM_NO_ERROR)
				CLEANUP("indexing failed", ret);

			if (type != s->type)
				CLEANUP("failed on type", -1);

			if (fabs(lattice_constant - 2) > tolerance)
				CLEANUP("failed on lattice constant", -1);

			double x = pdata[i][1][0];
			double y = pdata[i][1][1];
			double z = pdata[i][1][2];
			double iad = sqrt(x*x + y*y + z*z);
			if (fabs(iad - interatomic_distance) > tolerance)
				CLEANUP("failed on interatomic distance", -1);

			num_tests++;
		}
	}

cleanup:
	printf("num tests completed: %d\n", num_tests);
	ptm_uninitialize_local(local_handle);
	return ret;
}

