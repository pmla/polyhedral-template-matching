#include <cmath>
#include <cstdint>
#include <cstdbool>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <cassert>
#include "ptm_polar.h"
#include "ptm_arb_handle.h"
#include "ptm_arb_register.h"


namespace ptm {

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))


static double evaluate(int num_fixed, uint8_t* permutation_P, uint8_t* permutation_Q, double G1, double G2, double* A, double (*P)[3], double (*Q)[3], double* rot)
{
	if (num_fixed == 1)
	{
		double x1 = P[permutation_P[0]][0];
		double y1 = P[permutation_P[0]][1];
		double z1 = P[permutation_P[0]][2];

		double x2 = Q[permutation_Q[0]][0];
		double y2 = Q[permutation_Q[0]][1];
		double z2 = Q[permutation_Q[0]][2];

		double norm1 = sqrt(x1*x1 + y1*y1 + z1*z1);
		double norm2 = sqrt(x2*x2 + y2*y2 + z2*z2);
		double d = norm1 - norm2;
		return d*d;
	}

	double E0 = (G1 + G2) / 2;
	double nrmsdsq, q[4];
	FastCalcRMSDAndRotation(A, E0, &nrmsdsq, q, rot);
	assert(nrmsdsq == nrmsdsq);
	return nrmsdsq;
}

static void matvec(double* A, double* x, double* b)
{
	b[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
	b[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
	b[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}



/*
		double rotated_point[3] = {0, 0, 0};
		matvec(rotation, P[permutation_P[i]], rotated_point);

		double x0 = rotated_point[0];
		double y0 = rotated_point[1];
		double z0 = rotated_point[2];

		double x1 = Q[permutation_Q[i]][0];
		double y1 = Q[permutation_Q[i]][1];
		double z1 = Q[permutation_Q[i]][2];

		double dx = x1 - x0;
		double dy = y1 - y0;
		double dz = z1 - z0;
		double dist = dx*dx + dy*dy + dz*dz;

		acc += dist;
	}

	return sqrt(acc / num);

*/

#if 0
static int select_branch(int num_points, double (*P)[3], double (*Q)[3], uint8_t* permutation_P, uint8_t* permutation_Q, int level, double* rot, bool* p_branch_on_Q)
{
	double rotated_points[PTM_ARB_MAX_POINTS][3];
	for (int i=0;i<num_points;i++)
		matvec(rot, P[permutation_P[i]], rotated_points[i]);

	double d[PTM_ARB_MAX_POINTS][PTM_ARB_MAX_POINTS];
	for (int j=0;j<num_points;j++)
	{
		for (int i=0;i<num_points;i++)
		{
			double x1 = rotated_points[j][0];
			double y1 = rotated_points[j][1];
			double z1 = rotated_points[j][2];

			double x2 = Q[permutation_Q[i]][0];
			double y2 = Q[permutation_Q[i]][1];
			double z2 = Q[permutation_Q[i]][2];

			double dx = x2 - x1;
			double dy = y2 - y1;
			double dz = z2 - z1;

			double dist = dx*dx + dy*dy + dz*dz;
			d[j][i] = dist;
		}
	}

	int bi = -1;
	double max = 0.0;
	/*for (int j=level;j<num_points;j++)
	{
		double m = INFINITY;
		for (int i=level;i<num_points;i++)
			m = MIN(m, d[j][i]);

		if (m > max)
		{
			max = m;
			bi = j;
		}
	}*/

	for (int j=level;j<num_points;j++)
	{
		double m = INFINITY;
		for (int i=level;i<num_points;i++)
			m = MIN(m, d[i][j]);

		if (m > max)
		{
			max = m;
			bi = j;// + num_points;
		}
	}

	//*p_branch_on_Q = bi >= num_points;
	return bi % num_points;
}
#endif

static int _register_points(int num_points, double (*P)[3], double (*Q)[3], double vmax, bool quick, double* relaxation, uint8_t* best_permutation_P, uint8_t* best_permutation_Q, double* p_best, int* p_num_nodes_explored)
{
	handle_t* handle = handle_init(num_points, (size_t)1 * 1024*1024*1024);
	if (handle == NULL)
		return -1;

	double v = 0;
	int level = num_points;
	int it = 0;
	int num_nodes = 0;
	double best = INFINITY;
	double G1 = 0.0, G2 = 0.0, A[9] = {0};

	uint8_t permutation_P[PTM_ARB_MAX_POINTS];
	uint8_t permutation_Q[PTM_ARB_MAX_POINTS];
	for (int i=0;i<num_points;i++)
	{
		permutation_P[i] = i;
		permutation_Q[i] = i;
	}

	bool branch_on_Q = true;
	double rot[9];

	//push initial node (empty permutation)
	int ret = handle_add(handle, 0, 0.0, permutation_P, permutation_Q, 0.0, 0.0, A, branch_on_Q, 0);
	if (ret != 0)
		goto cleanup;

	//push best node (identity permutation)
	//int level = num_points;
	full_innerproduct(A, level, P, Q, permutation_P, permutation_Q, &G1, &G2);
	v = evaluate(level, permutation_P, permutation_Q, G1, G2, A, P, Q, rot) + relaxation[level];
	if (quick)
		v /= level;


	ret = handle_add(handle, level, v, permutation_P, permutation_Q, G1, G2, A, branch_on_Q, level);
	if (ret != 0)
		goto cleanup;

	best = v;

	if (quick)
	{
		memcpy(best_permutation_P, permutation_P, num_points * sizeof(uint8_t));
		memcpy(best_permutation_Q, permutation_Q, num_points * sizeof(uint8_t));
	}

	num_nodes = 2;
	while (pqueue_size(handle->pq) > 0)
	{
		int64_t index = 0;
		int ret = pqueue_pop(handle->pq, &v, &index);
		if (ret != 0)
			goto cleanup;

		node_t node;
		memcpy(&node, &handle->data[index], sizeof(node_t));
		if (v > best || v > vmax)
			continue;

		if (node.level < num_points)
		{
			memcpy(permutation_P, &handle->pscratch[node.permutation_index], num_points * sizeof(uint8_t));
			memcpy(permutation_Q, &handle->qscratch[node.permutation_index], num_points * sizeof(uint8_t));

			int num_swaps = num_points - node.level;
			for (int i=0;i<num_swaps;i++, it++)
			{
				if (node.branch_on_Q)
				{
					uint8_t temp = permutation_P[node.level];
					permutation_P[node.level] = permutation_P[node.level + i];
					permutation_P[node.level + i] = temp;
				}
				else
				{
					uint8_t temp = permutation_Q[node.level];
					permutation_Q[node.level] = permutation_Q[node.level + i];
					permutation_Q[node.level + i] = temp;
				}

				level = node.level + 1;
				G1 = node.G1;
				G2 = node.G2;
				memcpy(A, node.A, 9 * sizeof(double));
				increment_innerproduct(A, node.level, P, Q, permutation_P, permutation_Q, &G1, &G2);

				num_nodes++;
				v = evaluate(level, permutation_P, permutation_Q, G1, G2, A, P, Q, rot) + relaxation[level];
				if (quick)
					v /= level;


				if (v < best && v < vmax)
				{
					if (level == num_points)
					{
						best = v;
						memcpy(best_permutation_P, permutation_P, num_points * sizeof(uint8_t));
						memcpy(best_permutation_Q, permutation_Q, num_points * sizeof(uint8_t));
					}
					else
					{
						int branch_index = level;
						branch_on_Q = true;
						if (!quick)
						{
							//branch_index = select_branch(num_points, P, Q, permutation_P, permutation_Q, level, rot, &branch_on_Q);
						}

						//printf("level: %d v: %f / %f\n", level, v, best);
						int ret = handle_add(handle, level, v, permutation_P, permutation_Q, G1, G2, A, branch_on_Q, branch_index);
						if (ret != 0)
							goto cleanup;
					}

				}

				if (node.branch_on_Q)
				{
					uint8_t temp = permutation_P[node.level];
					permutation_P[node.level] = permutation_P[node.level + i];
					permutation_P[node.level + i] = temp;
				}
				else
				{
					uint8_t temp = permutation_Q[node.level];
					permutation_Q[node.level] = permutation_Q[node.level + i];
					permutation_Q[node.level + i] = temp;
				}
			}
		}
	}

cleanup:
	*p_best = best;
	*p_num_nodes_explored = num_nodes;
	handle_uninit(&handle);
	return ret;
}

typedef struct
{
	double x;
	int index;
} sorthelper_t;

static int sorthelper_compare(const void* va, const void* vb)
{
	sorthelper_t* a = (sorthelper_t*)va;
	sorthelper_t* b = (sorthelper_t*)vb;

	return a->x < b->x;
}

static int double_compare(const void* va, const void* vb)
{
	double a = *(double*)va;
	double b = *(double*)vb;
	return a > b;
}

static int sort_points_by_uniqueness(int num_points, double (*P)[3], double (*sorted_points)[3])
{
	double distances[PTM_ARB_MAX_POINTS][PTM_ARB_MAX_POINTS];
	for (int i=0;i<num_points;i++)
	{
		for (int j=0;j<num_points;j++)
		{
			double x1 = P[i][0];
			double y1 = P[i][1];
			double z1 = P[i][2];

			double x2 = P[j][0];
			double y2 = P[j][1];
			double z2 = P[j][2];

			double x = x1 - x2;
			double y = y1 - y2;
			double z = z1 - z2;
			distances[i][j] = sqrt(x*x + y*y + z*z);
		}
	}

	for (int i=0;i<num_points;i++)
		qsort(&distances[i], num_points, sizeof(double), double_compare);

	sorthelper_t furthest[PTM_ARB_MAX_POINTS];
	for (int i=0;i<num_points;i++)
	{
		double min_dist = INFINITY;
		for (int j=0;j<num_points;j++)
		{
			if (i == j) continue;

			double acc = 0.0;
			for (int k=0;k<num_points;k++)
			{
				double x = distances[i][k] - distances[j][k];
				acc += x*x;
			}

			min_dist = MIN(min_dist, acc);
		}

		furthest[i].index = i;
		furthest[i].x = min_dist;
	}

	qsort(furthest, num_points, sizeof(sorthelper_t), sorthelper_compare);

	for (int i=0;i<num_points;i++)
		memcpy(&sorted_points[i], &P[furthest[i].index], 3 * sizeof(double));
//return 0;

	int indices[PTM_ARB_MAX_POINTS];
	indices[0] = furthest[0].index;

	int bi = -1;
	double max = 0;
	for (int i=0;i<num_points;i++)
	{
		double x = P[i][0];
		double y = P[i][1];
		double z = P[i][2];
		double norm = x*x + y*y + z*z;
		if (norm > max)
		{
			max = norm;
			bi = i;
		}
	}

	indices[0] = bi;


	max = 0;
	/*
	int a = -1, b = -1;
	for (int i=0;i<num_points;i++)
	{
		for (int j=0;j<num_points;j++)
		{
			if (distances[i][j] > max)
			{
				max = distances[i][j];
				a = i;
				b = j;
			}
		}
	}

	//indices[0] = a;
	//indices[1] = b;
	*/

	int num_hit = 1;
	for (int i=num_hit;i<num_points;i++)
	{
		int bii = -1;
		double max = 0.0;
		for (int j=0;j<num_points;j++)
		{
			bool skip = false;
			for (int k=0;k<num_hit;k++)
				if (j == indices[k])
					skip = true;

			if (skip)
				continue;

			double min = INFINITY;
			for (int k=0;k<num_hit;k++)
			{
				double x1 = P[j][0];
				double y1 = P[j][1];
				double z1 = P[j][2];

				double x2 = P[indices[k]][0];
				double y2 = P[indices[k]][1];
				double z2 = P[indices[k]][2];

				double x = x1 - x2;
				double y = y1 - y2;
				double z = z1 - z2;
				double dist = sqrt(x*x + y*y + z*z);
				min = MIN(min, dist);
			}

			if (min > max)
			{
				bii = j;
				max = min;
			}
		}

		assert(bi != -1);
		indices[num_hit++] = bii;
	}

	assert(num_hit == num_points);
	for (int i=0;i<num_points;i++)
		memcpy(&sorted_points[i], &P[indices[i]], 3 * sizeof(double));

	return 0;
}

static double check_rmsd(int num, double (*P)[3], double (*Q)[3], uint8_t* permutation_P, uint8_t* permutation_Q, double* rotation)
{
	double acc = 0;
	for (int i=0;i<num;i++)
	{
		double rotated_point[3] = {0, 0, 0};
		matvec(rotation, P[permutation_P[i]], rotated_point);

		double x0 = rotated_point[0];
		double y0 = rotated_point[1];
		double z0 = rotated_point[2];

		double x1 = Q[permutation_Q[i]][0];
		double y1 = Q[permutation_Q[i]][1];
		double z1 = Q[permutation_Q[i]][2];

		double dx = x1 - x0;
		double dy = y1 - y0;
		double dz = z1 - z0;
		double dist = dx*dx + dy*dy + dz*dz;

		acc += dist;
	}

	return sqrt(acc / num);
}

int register_points(int num_points, double (*P)[3], double (*_Q)[3], double max_rmsd, uint8_t* best_permutation, double* p_rmsd, int* p_num_nodes_explored)
{
	assert(num_points <= PTM_ARB_MAX_POINTS);

	double Q[PTM_ARB_MAX_POINTS][3];
	if (1)
	{
		sort_points_by_uniqueness(num_points, _Q, Q);
	}
	else
	{
		memcpy(Q, _Q, 3 * num_points * sizeof(double));
	}

//for (int i=0;i<num_points;i++)
//	printf("%f %f %f\n", P[i][0] - Q[i][0], P[i][1] - Q[i][1], P[i][2] - Q[i][2]);

	double P_norms[PTM_ARB_MAX_POINTS];
	double Q_norms[PTM_ARB_MAX_POINTS];
	double relaxation[PTM_ARB_MAX_POINTS + 1];
	memset(relaxation, 0, (num_points + 1) * sizeof(double));
	for (int i=0;i<num_points;i++)
	{
		double x = P[i][0];
		double y = P[i][1];
		double z = P[i][2];
		P_norms[i] = sqrt(x*x + y*y + z*z);
	}

	for (int i=0;i<num_points;i++)
	{
		double x = Q[i][0];
		double y = Q[i][1];
		double z = Q[i][2];
		Q_norms[i] = sqrt(x*x + y*y + z*z);
	}

	for (int i=0;i<num_points;i++)
	{
		double smin = INFINITY;
		for (int j=0;j<num_points;j++)
		{
			double d = P_norms[j] - Q_norms[i];
			smin = MIN(smin, d*d);
		}

		relaxation[i] = smin;
	}

	for (int i=0;i<num_points;i++)
		for (int j=i+1;j<num_points;j++)
			relaxation[i] += relaxation[j];


	uint8_t best_permutation_P[PTM_ARB_MAX_POINTS];
	uint8_t best_permutation_Q[PTM_ARB_MAX_POINTS];

	int ret = 0, num_nodes0 = 0, num_nodes1 = 0;
	double vmax = max_rmsd*max_rmsd, best0 = INFINITY;
	ret = _register_points(num_points, P, Q, vmax, true, relaxation, best_permutation_P, best_permutation_Q, &best0, &num_nodes0);
	if (ret != 0)
		return ret;

	//printf("\n        %f\tnum nodes: %d\n", sqrt(best0), num_nodes0);
//printf("vmax: %f %f\n", vmax * num_points, best0 * num_points);

	vmax = MIN(vmax * num_points, best0 * num_points);

	double best1 = INFINITY;
	ret = _register_points(num_points, P, Q, vmax, false, relaxation, best_permutation_P, best_permutation_Q, &best1, &num_nodes1);
	if (ret != 0)
		return ret;

	best1 /= num_points;
	double best = MIN(best0, best1);

	//printf("bests: %f %f %f\n", vmax, best0, best1);
	//printf("\t%d %d\n", num_nodes0, num_nodes1);


	{
		/*printf("best permutation P:  ");
		for (int i=0;i<num_points;i++)
			printf("%d ", best_permutation_P[i]);
		printf("\n");
		printf("best permutation Q:  ");
		for (int i=0;i<num_points;i++)
			printf("%d ", best_permutation_Q[i]);
		printf("\n");*/

		double best_rotation[9];
		double G1 = 0.0, G2 = 0.0, A[9] = {0};
		int level = num_points;
		full_innerproduct(A, level, P, Q, best_permutation_P, best_permutation_Q, &G1, &G2);
		double v = evaluate(level, best_permutation_P, best_permutation_Q, G1, G2, A, P, Q, best_rotation);
		//printf("%f %f %f %f %f %f %f %f %f\n", best_rotation[0], best_rotation[1], best_rotation[2], best_rotation[3], best_rotation[4], best_rotation[5], best_rotation[6], best_rotation[7], best_rotation[8]);

		double _rmsd = check_rmsd(num_points, P, Q, best_permutation_P, best_permutation_Q, best_rotation);
		//printf("rmsd: %f %f\n", sqrt(best), _rmsd);
		//printf("v: %f\n", v);
		assert(fabs(_rmsd - sqrt(best)) < 1E-5);
	}

	*p_rmsd = sqrt(best);
	*p_num_nodes_explored = num_nodes0 + num_nodes1;
	return ret;
}
//todo: proper permutation

}

