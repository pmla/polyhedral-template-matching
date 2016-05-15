/*
	add icosahedral fundamental zone rotations

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

	cubic materials
		orientations + fundamental zones
		strain rotation must equal to rmsd rotation if no anisotropic distortion, when mapped into fundamental zone

	alloys
		disordered + ordered

	deformation gradients

	strains + polar decomposition rotations

	check det(U) > 0
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include "index_ptm.h"
#include "normalize_vertices.h"
#include "qcprot/quat.h"


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#define RADIANS(x) (2.0 * M_PI * (x) / 360.0)
#define DEGREES(x) (360 * (x) / (2.0 * M_PI))


typedef struct
{
	int type;
	int check;
	int num_points;
	const double (*points)[3];
} structdata_t;

structdata_t structdata[5] =  {	{ .type = PTM_MATCH_SC,  .check = PTM_CHECK_SC,  .num_points =  7, .points = ptm_template_sc },
				{ .type = PTM_MATCH_FCC, .check = PTM_CHECK_FCC, .num_points = 13, .points = ptm_template_fcc},
				{ .type = PTM_MATCH_HCP, .check = PTM_CHECK_HCP, .num_points = 13, .points = ptm_template_hcp},
				{ .type = PTM_MATCH_ICO, .check = PTM_CHECK_ICO, .num_points = 13, .points = ptm_template_ico},
				{ .type = PTM_MATCH_BCC, .check = PTM_CHECK_BCC, .num_points = 15, .points = ptm_template_bcc}};

typedef struct
{
	int32_t type;
	int32_t numbers[15];
} alloytest_t;

alloytest_t sc_alloy_tests[] = {

	{.type = PTM_ALLOY_NONE,   .numbers = {-1}},	//no test
};

alloytest_t fcc_alloy_tests[] = {

	{.type = PTM_ALLOY_NONE,   .numbers = {-1}},	//no test

	{.type = PTM_ALLOY_NONE,   .numbers = {0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},	//pure -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 4, 4, 4}},	//pure -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 7}},	//L10 -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {3, 0, 3, 0, 0, 3, 3, 3, 3, 0, 0, 0, 0}},	//L10 -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {3, 0, 1, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3}},	//L10 -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {0, 3, 1, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0}},	//L12_CU -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {5, 5, 3, 5, 5, 3, 3, 3, 3, 5, 5, 5, 5}},	//L12_CU -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {9, 9, 2, 9, 9, 9, 9, 9, 9, 3, 3, 3, 3}},	//L12_CU -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {4, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},	//L12_AU -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}},	//L12_AU -defect

	{.type = PTM_ALLOY_PURE,   .numbers = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
	{.type = PTM_ALLOY_PURE,   .numbers = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},
	{.type = PTM_ALLOY_L10,    .numbers = {3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0}},
	{.type = PTM_ALLOY_L10,    .numbers = {3, 0, 0, 0, 0, 3, 3, 3, 3, 0, 0, 0, 0}},
	{.type = PTM_ALLOY_L10,    .numbers = {3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3}},
	{.type = PTM_ALLOY_L12_CU, .numbers = {0, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0}},
	{.type = PTM_ALLOY_L12_CU, .numbers = {5, 5, 5, 5, 5, 3, 3, 3, 3, 5, 5, 5, 5}},
	{.type = PTM_ALLOY_L12_CU, .numbers = {9, 9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 3, 3}},
	{.type = PTM_ALLOY_L12_AU, .numbers = {4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
	{.type = PTM_ALLOY_L12_AU, .numbers = {1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}},
};

alloytest_t hcp_alloy_tests[] = {

	{.type = PTM_ALLOY_NONE,   .numbers = {-1}},	//no test
};

alloytest_t ico_alloy_tests[] = {

	{.type = PTM_ALLOY_NONE,   .numbers = {-1}},	//no test
};

alloytest_t bcc_alloy_tests[] = {

	{.type = PTM_ALLOY_NONE,   .numbers = {-1}},	//no test

	{.type = PTM_ALLOY_NONE,   .numbers = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}},	//pure -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},	//pure -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {4, 0, 0, 0, 0, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4}},	//B2 -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {1, 2, 2, 2, 2, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1}},	//B2 -defect

	{.type = PTM_ALLOY_PURE,   .numbers = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
	{.type = PTM_ALLOY_PURE,   .numbers = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},
	{.type = PTM_ALLOY_B2,     .numbers = {4, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4}},
	{.type = PTM_ALLOY_B2,     .numbers = {1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1}},
};

typedef struct
{
	double pre[9];
	double post[9];
	bool fundamental;
	bool strain;
} quattest_t;

quattest_t qtest[] = {	{	.pre  = {1.00, 0.00, -0.00, 0.00},
				.post = {1.00, 0.00, -0.00, 0.00},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.00, 1.00, -0.00, 0.00},
				.post = {1.00, 0.00, -0.00, 0.00},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.00, 0.00, -1.00, 0.00},
				.post = {1.00, 0.00, -0.00, 0.00},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.00, 0.00, -0.00, 1.00},
				.post = {1.00, 0.00, -0.00, 0.00},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.987070, 0.020780, -0.031171, 0.155853},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.683270, 0.712658, 0.088164, 0.132246},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.720005, -0.095511, 0.675923, 0.124899},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.587759, -0.007347, -0.036735, 0.808168},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.808168, 0.036735, -0.007347, -0.587759},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.675923, 0.124899, -0.720005, 0.095511},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.712658, -0.683270, -0.132246, 0.088164},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.420803, 0.410413, 0.545486, 0.597437},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.576656, 0.441584, 0.566266, -0.389633},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.389633, 0.566266, -0.441584, 0.576656},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.545486, 0.597437, -0.420803, -0.410413},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.441584, -0.576656, 0.389633, 0.566266},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.597437, -0.545486, 0.410413, -0.420803},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.410413, -0.420803, -0.597437, 0.545486},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.566266, -0.389633, -0.576656, -0.441584},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.020780, -0.987070, -0.155853, -0.031171},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.007347, 0.587759, 0.808168, 0.036735},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.124899, -0.675923, -0.095511, -0.720005},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.095511, 0.720005, 0.124899, -0.675923},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.036735, -0.808168, 0.587759, -0.007347},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.031171, -0.155853, 0.987070, 0.020780},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.088164, 0.132246, -0.683270, -0.712658},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.132246, -0.088164, 0.712658, -0.683270},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {0.155853, 0.031171, 0.020780, -0.987070},
				.post = {0.987070, 0.020780, -0.031171, 0.155853},
				.fundamental = true,
				.strain = false				},

			{	.pre  = {1.1, 0.0, 0.0, -0.0, 1.0, -0.0, -0.0, 0.0, 1.0},
				.post = {},
				.fundamental = false,
				.strain = true				},

			{	.pre  = {1.05, 0.0, 0.0, -0.05, 0.95, -0.0, -0.02, 0.06, 0.94},
				.post = {},
				.fundamental = false,
				.strain = true				},

			};

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
	double transformed_template[15][3];
	for (int i=0;i<num;i++)
	{
		double row[3] = {0, 0, 0};
		matvec(A, (double*)template_points[i], row);
		memcpy(transformed_template[i], row, 3 * sizeof(double));
	}

	//translate and scale input points
	double points[15][3];
	subtract_barycentre(num, input_points[0], points);
	for (int i=0;i<num;i++)
		for (int j=0;j<3;j++)
			points[i][j] *= scale;

	double acc = 0;
	for (int i=0;i<num;i++)
	{
		double x0 = points[i][0];
		double y0 = points[i][1];
		double z0 = points[i][2];

		int bi = -1;
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
			if (dist < min_dist)
				bi = j;
			min_dist = MIN(min_dist, dist);
		}

		acc += min_dist;
	}

	return sqrt(fabs(acc / num));
}

static double mapped_neighbour_rmsd(int num, double scale, double* A, double (*input_points)[3], const double (*template_points)[3], int8_t* mapping)
{
	//transform template
	double transformed_template[15][3];
	for (int i=0;i<num;i++)
	{
		double row[3] = {0, 0, 0};
		matvec(A, (double*)template_points[i], row);
		memcpy(transformed_template[i], row, 3 * sizeof(double));
	}

	//translate and scale input points
	double points[15][3];
	subtract_barycentre(num, input_points[0], points);
	for (int i=0;i<num;i++)
		for (int j=0;j<3;j++)
			points[i][j] *= scale;

	/*printf("\n");
	for (int i=0;i<num;i++)
		printf("%d ", mapping[i]);
	printf("\n");

	for (int i=0;i<num;i++)
		printf("%f %f %f\n", points[mapping[i]][0], points[mapping[i]][1], points[mapping[i]][2]);
	printf("\n");
	for (int i=0;i<num;i++)
		printf("%f %f %f\n", transformed_template[i][0], transformed_template[i][1], transformed_template[i][2]);
	printf("\n");*/

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

		//printf("dist: %f\n", dist);
		acc += dist;
	}

	return sqrt(fabs(acc / num));
}

uint64_t run_tests()
{
	int ret = 0;
	const double tolerance = 1E-5;
	double identity_matrix[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	int num_structures = sizeof(structdata) / sizeof(structdata_t);
	int num_quat_tests = sizeof(qtest) / sizeof(quattest_t);

	int num_alloy_tests[] = {	sizeof(sc_alloy_tests) / sizeof(alloytest_t),
					sizeof(fcc_alloy_tests) / sizeof(alloytest_t),
					sizeof(hcp_alloy_tests) / sizeof(alloytest_t),
					sizeof(ico_alloy_tests) / sizeof(alloytest_t),
					sizeof(bcc_alloy_tests) / sizeof(alloytest_t)	};

	alloytest_t* alloy_test[] = {	sc_alloy_tests,
					fcc_alloy_tests,
					hcp_alloy_tests,
					ico_alloy_tests,
					bcc_alloy_tests	};
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
int it = 4;
	{
		structdata_t* s = &structdata[it];
		double points[15][3];
		double rotated_points[15][3];

		memcpy(points, s->points, 3 * sizeof(double) * s->num_points);

		for (int k=0;k<24;k++)
		{
			double U[9];
			extern double generator_cubic[24][4];
			double qtemp[4];
			memcpy(qtemp, generator_cubic[k], 4 * sizeof(double));
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

	for (int iq=0;iq<num_quat_tests;iq++)
	{
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

		for (int it = 0;it<num_structures;it++)
		{
			structdata_t* s = &structdata[it];
			double points[15][3];
			memcpy(points, s->points, 3 * sizeof(double) * s->num_points);

			for (int i=0;i<s->num_points;i++)
			{
				double row[3] = {0, 0, 0};
				matvec(rot, points[i], row);
				memcpy(points[i], row, 3 * sizeof(double));
			}

			double rescale = 0.5;
			double offset[3] = {12.0, 45.6, 789.10};
			for (int i = 0;i<s->num_points;i++)
				for (int j = 0;j<3;j++)
					points[i][j] = points[i][j] * rescale + offset[j];

			int tocheck = 0;
			for (int i = 0;i<it+1;i++)
				tocheck |= structdata[i].check;

			for (int ia=0;ia<num_alloy_tests[it];ia++)
			{
				int32_t* numbers = alloy_test[it][ia].numbers;
				if (numbers[0] == -1)
					numbers = NULL;

				for (int itop=0;itop<=1;itop++)
				{
					int8_t mapping[15];
					bool topological = itop == 1;
					int32_t type, alloy_type;
					double scale, rmsd;
					double q[4], F[9], F_res[3], U[9], P[9];
					ret = ptm_index(local_handle, s->num_points, points[0], numbers, tocheck, topological, &type, &alloy_type, &scale, &rmsd, q, F, F_res, U, P, mapping);
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
#endif

					//check type
					if (type != s->type)
						CLEANUP("failed on type", -1);

					//check alloy type
					if (alloy_type != alloy_test[it][ia].type)
						CLEANUP("failed on alloy type", -1);

					if (matrix_determinant(U) <= 0)
						CLEANUP("failed on U-matrix right-handedness test", -1);

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
						if (!qtest[iq].fundamental || (s->type == PTM_MATCH_SC || s->type == PTM_MATCH_FCC || s->type == PTM_MATCH_BCC))
							if (quat_misorientation(q, qpost) > tolerance)
								CLEANUP("failed on disorientation", -1);

						//check deformation gradient disorientation
						if (s->type == PTM_MATCH_SC || s->type == PTM_MATCH_FCC || s->type == PTM_MATCH_BCC)
						{
							double qu[4];
							rotation_matrix_to_quaternion(U, qu);
							//printf("%f %f %f %f\n", qu[0], qu[1], qu[2], qu[3]);
							//printf("%f %f %f %f\n", qpost[0], qpost[1], qpost[2], qpost[3]);
							if (quat_misorientation(qu, qpost) > tolerance)
								CLEANUP("failed on deformation gradient disorientation", -1);
						}

						quaternion_to_rotation_matrix(q, A);
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

#ifdef DEBUG
#endif

cleanup:
	printf("num tests completed: %d\n", num_tests);
	ptm_uninitialize_local(local_handle);
	return ret;
}

