#include <cstdint>
#include <cstdbool>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "ptm_constants.h"
#include "ptm_quat.h"
#include "ptm_normalize_vertices.h"
#include "ptm_arb_register_dfs.h"
#include "ptm_arb_treefilter.h"


static double frand(double fmin, double fmax)
{
	double f = (double)rand() / RAND_MAX;
	return fmin + f * (fmax - fmin);
}

static void random_quaternion(double* q)
{
	double x0 = frand(0.0, 1.0);
	double x1 = frand(0.0, 2 * M_PI);
	double x2 = frand(0.0, 2 * M_PI);

	double r1 = sqrt(1 - x0);
	double r2 = sqrt(x0);

	double s1 = sin(x1);
	double c1 = cos(x1);
	double s2 = sin(x2);
	double c2 = cos(x2);

	q[0] = s1*r1;
	q[1] = c1*r1;
	q[2] = s2*r2;
	q[3] = c2*r2;
}

static void random_rotation_matrix(double* u)
{
	double q[4];
	random_quaternion(q);
	ptm::quaternion_to_rotation_matrix(q, u);
}

static void randomly_permute_points(int num_points, double (*points)[3], double (*permuted)[3], bool permute)
{
	int temp[num_points];
	int indices[num_points];

	for (int i=0;i<num_points;i++)
		temp[i] = i;

	int num_left = num_points;
	for (int i=0;i<num_points;i++)
	{
		int t = rand() % num_left;
		if (i == 0)
			t = 0;	//skip first point

		indices[i] = temp[t];
		temp[t] = temp[--num_left];
	}

	for (int i=0;i<num_points;i++)
		for (int j=0;j<3;j++)
			if (permute)
				permuted[i][j] = points[indices[i]][j];
			else
				permuted[i][j] = points[i][j];
}

static void matvec(double* A, double* x, double* b)
{
	b[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
	b[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
	b[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}

static void random_point_rotation(int num_points, double (*points)[3])
{
	double U[9];
	random_rotation_matrix(U);

	for (int i=0;i<num_points;i++)
	{
		double temp[3];
		matvec(U, points[i], temp);
		memcpy(points[i], temp, 3 * sizeof(double));
	}
}

int main()
{
	uint64_t seed = 0;
	srand(seed);

	int num_points = 23;
	assert(num_points <= 255);

	double P[num_points][3];
	double Q[num_points][3];

	const double k = 0.05;
	const double max_rmsd = 0.15;//9990.25;


	int total_nodes_astar = 0;
	double rmsd_sum_astar = 0.0;
	int total_nodes_dfs = 0;
	double rmsd_sum_dfs = 0.0;
	int num_it = 10000;
int num_bad = 0;
double worst = 0;
	for (int it=0;it<num_it;it++)
	{
		for (int i=0;i<num_points;i++)
			for (int j=0;j<3;j++)
				P[i][j] = frand(-1, 1);
		memcpy(P, ptm_template_fluorite_f, num_points * 3 * sizeof(double));

		for (int i=0;i<num_points;i++)
			for (int j=0;j<3;j++)
				Q[i][j] = frand(-1, 1);
		Q[0][0] = Q[0][1] = Q[0][2] = 0;

		randomly_permute_points(num_points, P, Q, true);
		random_point_rotation(num_points, Q);

		double m = 0.0;
		for (int i=0;i<num_points;i++)
		{
			for (int j=0;j<3;j++)
			{
				//P[i][j] = frand(-1, 1);

				double p[3] = {P[i][0], P[i][1], P[i][2]};
				Q[i][j] = P[i][j] + frand(-k, k);
				double dx = p[0] - P[i][0];
				double dz = p[1] - P[i][1];
				double dy = p[2] - P[i][2];
				m += dx*dx + dy*dy + dz*dz;
			}
		}
		//printf("m: %f\n", m);
	
		//ptm::normalize_vertices(num_points, P, P)
		ptm::normalize_vertices(num_points, Q, Q);	//todo: save scale

		double rmsd_astar = 0.0;
		int num_nodes_astar = 0;
		uint8_t best_permutation[num_points];
		int ret = 0;//ptm::register_points(num_points, P, Q, max_rmsd, best_permutation, &rmsd_astar, &num_nodes_astar);
		if (ret != 0)
		{
			printf("registration failed\n");
			exit(2);
		}

		if (rmsd_astar > worst)
			num_bad++;
		worst = std::max(worst, rmsd_astar);


		int num_nodes_dfs = 0;
		double rmsd_dfs = INFINITY;
		int match_found = ptm::register_points_dfs(num_points, P, Q, max_rmsd, &ptm::filter_fluorite_f,
						best_permutation, &rmsd_dfs, &num_nodes_dfs);
		if (match_found)
			rmsd_sum_dfs += rmsd_dfs;
		total_nodes_dfs += num_nodes_dfs;

		printf("%d rmsd: %f (%f)\tnum nodes: %d (%d)\n", it, rmsd_astar, rmsd_dfs, num_nodes_astar, num_nodes_dfs);
		total_nodes_astar += num_nodes_astar;
		rmsd_sum_astar += rmsd_astar;

		//if (fabs(rmsd_astar - rmsd_dfs) > 1E-4)
		//	printf("mismatch: %f %f\n", rmsd_astar, rmsd_dfs);
	}

printf("num_bad: %d %f\n", num_bad, worst);

	printf("total nodes astar: %d\n", total_nodes_astar);
	printf("average astar: %f\n", (double)total_nodes_astar / num_it);
	printf("rmsd average: %f\n", rmsd_sum_astar / num_it);

	printf("\n");
	printf("total nodes dfs: %d\n", total_nodes_dfs);
	printf("average dfs: %f\n", (double)total_nodes_dfs / num_it);
	printf("rmsd average dfs: %f\n", rmsd_sum_dfs / num_it);

	printf("difference: %f\n", rmsd_sum_dfs - rmsd_sum_astar);
	return 0;
}

