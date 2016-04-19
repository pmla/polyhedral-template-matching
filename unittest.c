#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include "index_PTM.h"
#include "quat.h"
#include "normalize_vertices.h"


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#define RADIANS(x) (2.0 * M_PI * (x) / 360.0)
#define DEGREES(x) (360 * (x) / (2.0 * M_PI))


extern const double points_sc[7][3];
extern const double points_fcc[13][3];
extern const double points_hcp[13][3];
extern const double points_ico[13][3];
extern const double points_bcc[15][3];

typedef struct
{
	int type;
	int check;
	int num_points;
	const double (*points)[3];
} structdata_t;

structdata_t structdata[5] =  {	{ .type = PTM_MATCH_SC,  .check = PTM_CHECK_SC,  .num_points =  7, .points = points_sc },
				{ .type = PTM_MATCH_FCC, .check = PTM_CHECK_FCC, .num_points = 13, .points = points_fcc},
				{ .type = PTM_MATCH_HCP, .check = PTM_CHECK_HCP, .num_points = 13, .points = points_hcp},
				{ .type = PTM_MATCH_ICO, .check = PTM_CHECK_ICO, .num_points = 13, .points = points_ico},
				{ .type = PTM_MATCH_BCC, .check = PTM_CHECK_BCC, .num_points = 15, .points = points_bcc}};


/*
	add icosahedral fundamental zone rotations

	structures:
		+ no match for each of them
		nonzero rmsd

	deformation gradients
		+ residuals

	strains + polar decomposition rotations

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
*/

typedef struct
{
	double pre[4];
	double post[4];
	bool fundamental;
} quattest_t;

typedef struct
{
	int32_t type;
	int32_t numbers[15];
} alloytest_t;

#define NUM_STRUCTURES 5
#define NUM_QUAT_TESTS 5
#define NUM_FCC_TESTS 20
#define NUM_BCC_TESTS 8

alloytest_t fcc_alloy_tests[NUM_FCC_TESTS] = {

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

alloytest_t bcc_alloy_tests[NUM_BCC_TESTS] = {

	{.type = PTM_ALLOY_NONE,   .numbers = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}},	//pure -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},	//pure -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {4, 0, 0, 0, 0, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4}},	//B2 -defect
	{.type = PTM_ALLOY_NONE,   .numbers = {1, 2, 2, 2, 2, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1}},	//B2 -defect

	{.type = PTM_ALLOY_PURE,   .numbers = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
	{.type = PTM_ALLOY_PURE,   .numbers = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}},
	{.type = PTM_ALLOY_B2,     .numbers = {4, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4}},
	{.type = PTM_ALLOY_B2,     .numbers = {1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1}},
};

quattest_t qtest[NUM_QUAT_TESTS] = {	{	.pre  = {1.00, 0.00, -0.00, 0.00},
						.post = {1.00, 0.00, -0.00, 0.00},
						.fundamental = true				},

					{	.pre  = {0.00, 1.00, -0.00, 0.00},
						.post = {1.00, 0.00, -0.00, 0.00},
						.fundamental = true				},

					{	.pre  = {0.00, 0.00, -1.00, 0.00},
						.post = {1.00, 0.00, -0.00, 0.00},
						.fundamental = true				},

					{	.pre  = {0.00, 0.00, -0.00, 1.00},
						.post = {1.00, 0.00, -0.00, 0.00},
						.fundamental = true				},

					{	.pre  = {0.95, 0.02, -0.03, 0.15},
						.post = {0.95, 0.02, -0.03, 0.15},
						.fundamental = true				},
					};


void matvec(double* A, double* x, double* b)
{
	b[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
	b[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
	b[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}

static double nearest_neighbour_rmsd(int num, double scale, double* q, double (*input_points)[3], const double (*template)[3])
{
	//rotate template
	double U[9];
	double rotated_template[15][3];
	quaternion_to_rotation_matrix(q, U);
	for (int i=0;i<num;i++)
	{
		double row[3] = {0, 0, 0};
		matvec(U, (double*)template[i], row);
		memcpy(rotated_template[i], row, 3 * sizeof(double));
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

		double min_dist = INFINITY;
		for (int j=0;j<num;j++)
		{
			double x1 = rotated_template[j][0];
			double y1 = rotated_template[j][1];
			double z1 = rotated_template[j][2];

			double dx = x1 - x0;
			double dy = y1 - y0;
			double dz = z1 - z0;
			double dist = dx*dx + dy*dy + dz*dz;
			min_dist = MIN(min_dist, dist);
		}

		acc += min_dist;
	}

	return sqrt(fabs(acc / num));
}

uint64_t run_tests()
{
	uint64_t result = 0;
	uint64_t bitmask = 1;
	const double tolerance = 1E-5;

	//rotation matrix => quaternion
	{
		double U[9] = {	 0,  0,  1,
				 0,  1,  0,
				-1,  0,  0  };
		double q[4];
		rotation_matrix_to_quaternion(U, q);
		double ans[4] = {1 / sqrt(2), 0, 1 / sqrt(2), 0};
		for (int i = 0;i<4;i++)
			assert(fabs(q[i] - ans[i]) < tolerance);
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
			assert(fabs(U[i] - ans[i]) < tolerance);
	}

	assert(NUM_QUAT_TESTS * NUM_STRUCTURES < 64);

	for (int iq=0;iq<NUM_QUAT_TESTS;iq++)
	{
		double qpre[4], qpost[4];
		memcpy(qpre, qtest[iq].pre, 4 * sizeof(double));
		memcpy(qpost, qtest[iq].post, 4 * sizeof(double));

		normalize_quaternion(qpre);
		normalize_quaternion(qpost);

		for (int it = 0;it<NUM_STRUCTURES;it++, bitmask <<= 1)
		{
			structdata_t* s = &structdata[it];
			double points[15][3];
			memcpy(points, s->points, 3 * sizeof(double) * s->num_points);

			double rot[9];
			quaternion_to_rotation_matrix(qpre, rot);

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
				tocheck |= s->check;

			alloytest_t* alloytest = NULL;
			int num_alloy_tests = 1;
			if (s->type == PTM_MATCH_FCC)
			{
				num_alloy_tests = NUM_FCC_TESTS;
				alloytest = &fcc_alloy_tests[0];
			}
			else if (s->type == PTM_MATCH_BCC)
			{
				num_alloy_tests = NUM_BCC_TESTS;
				alloytest = &bcc_alloy_tests[0];
			}

			for (int ia=0;ia<num_alloy_tests;ia++)
			{
				for (int itop=0;itop<=1;itop++)
				{
					bool topological = itop == 1;
					int32_t type, alloy_type;
					double scale, rmsd;
					double q[4], F[9], F_res[3], U[9], P[9];
					if (num_alloy_tests == 1)
						index_PTM(s->num_points, points[0], NULL, tocheck, topological, &type, &alloy_type, &scale, &rmsd, q, F, F_res, U, P);
					else
						index_PTM(s->num_points, points[0], alloytest[ia].numbers, tocheck, topological, &type, &alloy_type, &scale, &rmsd, q, F, F_res, U, P);

#ifdef DEBUG
					printf("type:\t\t%d\t(should be: %d)\n", type, s->type);
					printf("alloy type:\t%d\n", alloy_type);
					printf("scale:\t\t%f\n", scale);
					printf("rmsd:\t\t%f\n", rmsd);
					printf("quat: \t\t%.4f %.4f %.4f %.4f\n", q[0], q[1], q[2], q[3]);
					printf("qpost:\t\t%.4f %.4f %.4f %.4f\n", qpost[0], qpost[1], qpost[2], qpost[3]);
					//printf("rot: %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", rot[0], rot[1], rot[2], rot[3], rot[4], rot[5], rot[6], rot[7], rot[8]);
					//printf("U:   %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7], U[8]);
					//printf("P:   %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8]);
					//printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", rot[0] - U[0], rot[1] - U[1], rot[2] - U[2], rot[3] - U[3], rot[4] - U[4], rot[5] - U[5], rot[6] - U[6], rot[7] - U[7], rot[8] - U[8]);
#endif
					if (type != s->type)
					{
						result |= bitmask;
#ifdef DEBUG
						printf("failed on type\n");
#endif
					}

					if (num_alloy_tests > 1)
					{
						if ((s->type == PTM_MATCH_FCC || s->type == PTM_MATCH_BCC) && alloy_type != alloytest[ia].type)
						{
							result |= bitmask;
#ifdef DEBUG
							printf("failed on alloy type\n");
							printf("ia: %d\n", ia);
#endif
						}
					}

					if (rmsd > tolerance)
					{
						result |= bitmask;
#ifdef DEBUG
						printf("failed on rmsd\n");
#endif
					}

					if (fabs(scale - 1 / rescale) > tolerance)
					{
						result |= bitmask;
#ifdef DEBUG
						printf("failed on scale\n");
#endif
					}

					for (int i=0;i<3;i++)
					{
						for (int j=0;j<3;j++)
						{
							if (i == j && fabs(P[i*3 + j] - 1) > tolerance)
							{
								result |= bitmask;
#ifdef DEBUG
								printf("failed on P matrix diagonal\n");
#endif
							}
							else if (i != j && fabs(P[i*3 + j]) > tolerance)
							{
								result |= bitmask;
#ifdef DEBUG
								printf("failed on P matrix off-diagonal\n");
#endif
							}
						}
					}

					if (!qtest[iq].fundamental || (s->type == PTM_MATCH_SC || s->type == PTM_MATCH_FCC || s->type == PTM_MATCH_BCC))
					{
						if (quat_misorientation(q, qpost) > tolerance)
						{
							result |= bitmask;
#ifdef DEBUG
							printf("failed on disorientation\n");
#endif
						}
					}

					if (s->type == PTM_MATCH_SC || s->type == PTM_MATCH_FCC || s->type == PTM_MATCH_BCC)
					{
						double qu[4];
						rotation_matrix_to_quaternion(U, qu);
						rotate_quaternion_into_cubic_fundamental_zone(qu);

						if (quat_misorientation(qu, qpost) > tolerance)
						{
							result |= bitmask;
#ifdef DEBUG
							printf("failed on deformation gradient disorientation\n");
#endif
						}
					}

					double rmsd_approx = nearest_neighbour_rmsd(s->num_points, scale, q, points, s->points);
					if (fabs(rmsd_approx) > tolerance)
					{
						result |= bitmask;
#ifdef DEBUG
						printf("failed on rmsd nearest neighbour\n");
#endif
					}

#ifdef DEBUG
					printf("result: %lu\n", result);
					printf("\n\n");
					if (result != 0)
						return result;
#endif
				}
			}
		}
	}

	return result;
}

