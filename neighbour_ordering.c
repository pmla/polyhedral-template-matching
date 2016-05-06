#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define TOLERANCE 1E-8

#define MAX_POINTS 19
#define MAXF (2*MAX_POINTS - 4)


static double norm_squared_4D(double* p)
{
	double x = p[0];
	double y = p[1];
	double z = p[2];
	double w = p[3];

	return x*x + y*y + z*z + w*w;
}

static double vector_distance_2_4D(double* a, double* b)
{
	double x = a[0] - b[0];
	double y = a[1] - b[1];
	double z = a[2] - b[2];
	double w = a[3] - b[3];

	return x*x + y*y + z*z + w*w;
}

static double norm_squared(double* p)
{
	double x = p[0];
	double y = p[1];
	double z = p[2];

	return x*x + y*y + z*z;
}

static double dot_product(const double* a, const double* b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void cross_product(double* a, double* b, double* c)
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

static double dethelper(double x2, double x3, double x4, double y2, double y3, double y4, double z2, double z3, double z4)
{
	return - x2*y3*z4 + x2*y4*z3 + x3*y2*z4 - x3*y4*z2 - x4*y2*z3 + x4*y3*z2;
}

static bool circumsphere_centre(double* b, double* c, double* d, double* centre)
{
	//formula adapted from http://mathworld.wolfram.com/Circumsphere.html


	//norm squared is stored in fourth element
	double sb = b[3];
	double sc = c[3];
	double sd = d[3];

	double det = dethelper(b[0], c[0], d[0], b[1], c[1], d[1], b[2], c[2], d[2]);
	if (fabs(det) < 1E-6)
		return false;	//degenerate facet

	double detx = dethelper(sb, sc, sd, b[1], c[1], d[1], b[2], c[2], d[2]);
	double dety = dethelper(sb, sc, sd, b[0], c[0], d[0], b[2], c[2], d[2]);
	double detz = dethelper(sb, sc, sd, b[0], c[0], d[0], b[1], c[1], d[1]);

	centre[0] = detx / (2 * det);
	centre[1] = -dety / (2 * det);
	centre[2] = detz / (2 * det);
	return true;
}

static void cross_product_4D(const double* r0, const double* r1, const double* r2, double* p)
{
	p[0] =	+r0[1] * (r1[2] * r2[3] - r1[3] * r2[2])
		-r0[2] * (r1[1] * r2[3] - r1[3] * r2[1])
		+r0[3] * (r1[1] * r2[2] - r1[2] * r2[1]);

	p[1] =	-r0[0] * (r1[2] * r2[3] - r1[3] * r2[2])
		+r0[2] * (r1[0] * r2[3] - r1[3] * r2[0])
		-r0[3] * (r1[0] * r2[2] - r1[2] * r2[0]);

	p[2] =	+r0[0] * (r1[1] * r2[3] - r1[3] * r2[1])
		-r0[1] * (r1[0] * r2[3] - r1[3] * r2[0])
		+r0[3] * (r1[0] * r2[1] - r1[1] * r2[0]);

	p[3] =	-r0[0] * (r1[1] * r2[2] - r1[2] * r2[1])
		+r0[1] * (r1[0] * r2[2] - r1[2] * r2[0])
		-r0[2] * (r1[0] * r2[1] - r1[1] * r2[0]);
}

static void calculate_plane_normal_4D(const double (*points)[4], int b, int c, int d, double* plane_normal)
{
	cross_product_4D(points[b], points[c], points[d], plane_normal);
	double norm = sqrt(norm_squared_4D(plane_normal));
	plane_normal[0] /= norm;
	plane_normal[1] /= norm;
	plane_normal[2] /= norm;
	plane_normal[3] /= norm;
}

static double point_plane_distance_4D(const double* w, const double* plane_cross)
{
	return	  plane_cross[0] * (- w[0])
		+ plane_cross[1] * (- w[1])
		+ plane_cross[2] * (- w[2])
		+ plane_cross[3] * (- w[3]);
}

static bool visible(const double* w, const double* plane_normal)
{
	return point_plane_distance_4D(w, plane_normal) > 0;
}

static void _add_facet(const double (*points)[4], int b, int c, int d, int8_t* facet, double* plane_normal, double* barycentre)
{
#ifdef DEBUG
	printf("adding: %d %d %d %d\n", 0, b, c, d);
	assert(b < c);
	assert(b < d);
	assert(c < d);

	double _centre[3];
	if (!circumsphere_centre((double*)points[b], (double*)points[c], (double*)points[d], _centre))
	{
			printf("coplanar facet added! %d %d %d %d\n", 0, b, c, d);
			exit(3);
	}
#endif

	facet[0] = b;
	facet[1] = c;
	facet[2] = d;

	calculate_plane_normal_4D(points, b, c, d, plane_normal);
	if (visible(barycentre, plane_normal))
	{
		plane_normal[0] = -plane_normal[0];
		plane_normal[1] = -plane_normal[1];
		plane_normal[2] = -plane_normal[2];
		plane_normal[3] = -plane_normal[3];
	}
}

static void calculate_plane_normal(const double (*points)[4], int a, int b, int c, double* plane_normal)
{
	double u[3] = {	points[b][0] - points[a][0],
			points[b][1] - points[a][1],
			points[b][2] - points[a][2]	};

	double v[3] = {	points[c][0] - points[a][0],
			points[c][1] - points[a][1],
			points[c][2] - points[a][2]	};

	cross_product(u, v, plane_normal);
	double norm = sqrt(norm_squared(plane_normal));
	plane_normal[0] /= norm;
	plane_normal[1] /= norm;
	plane_normal[2] /= norm;
}

static double point_plane_distance(const double* w, const double* plane_point, const double* plane_cross)
{
	return	  plane_cross[0] * (plane_point[0] - w[0])
		+ plane_cross[1] * (plane_point[1] - w[1])
		+ plane_cross[2] * (plane_point[2] - w[2]);
}

static bool find_third_point(int num_points, const double (*points)[4], int a, int b, int* p_c)
{
	const double* x1 = points[a];
	const double* x2 = points[b];

	double x2x1[3] = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
	double ns_x2x1 = norm_squared(x2x1);

	int bi = -1;
	double max_dist = 0.0;
	for (int i = 0;i<num_points;i++)
	{
		if (i == a || i == b)
			continue;

		const double* x0 = points[i];

		double x1x0[3] = {x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2]};
		double dot = dot_product(x1x0, x2x1);
		double dist = (norm_squared(x1x0) * ns_x2x1 - dot*dot) / ns_x2x1;

		if (dist > max_dist)
		{
			max_dist = dist;
			bi = i;
		}
	}

	*p_c = bi;
	return max_dist > TOLERANCE;
}

static bool find_fourth_point(int num_points, const double (*points)[4], int a, int b, int c, int* p_d)
{
	double plane_normal[3];
	calculate_plane_normal(points, a, b, c, plane_normal);


	int bi = -1;
	double max_dist = 0.0;
	for (int i = 0;i<num_points;i++)
	{
		if (i == a || i == b || i == c)
			continue;

		const double* x0 = points[i];
		double dist = fabs(point_plane_distance(x0, points[a], plane_normal));
		if (dist > max_dist)
		{
			max_dist = dist;
			bi = i;
		}
	}

	*p_d = bi;
	return max_dist > TOLERANCE;
}

static bool find_fifth_point(int num_points, const double (*points)[4], int a, int b, int c, int d, int* p_e)
{
	double plane_normal[4];
	calculate_plane_normal_4D(points, b, c, d, plane_normal);


	int bi = -1;
	double max_dist = 0.0;
	for (int i = 0;i<num_points;i++)
	{
		if (i == a || i == b || i == c || i == d)
			continue;

		double _centre[3];
		if (!circumsphere_centre((double*)points[b], (double*)points[c], (double*)points[i], _centre)) continue;
		if (!circumsphere_centre((double*)points[b], (double*)points[d], (double*)points[i], _centre)) continue;
		if (!circumsphere_centre((double*)points[c], (double*)points[d], (double*)points[i], _centre)) continue;

		const double* x0 = points[i];
		double dist = fabs(point_plane_distance_4D(x0, plane_normal));
		if (dist > max_dist)
		{
			max_dist = dist;
			bi = i;
		}
	}

	*p_e = bi;
	return max_dist > TOLERANCE;
}

static int initial_simplex(int num_points, const double (*points)[4], int* initial_vertices)
{
	int a = 0;
	int bi = -1;
	double max_dist = 0.0;
	for (int i=0;i<num_points;i++)
	{
		double dist = vector_distance_2_4D((double*)points[a], (double*)points[i]);
		if (dist > max_dist)
		{
			bi = i;
			max_dist = dist;
		}
	}

	int b = bi, c = -1, d = -1;
	if (!find_third_point(num_points, points, a, b, &c))
		return -2;

	if (!find_fourth_point(num_points, points, a, b, c, &d))
		return -3;

	if (b > c)
	{
		int temp = b;
		b = c;
		c = temp;
	}
	if (c > d)
	{
		int temp = c;
		c = d;
		d = temp;
	}
	if (b > c)
	{
		int temp = b;
		b = c;
		c = temp;
	}

	initial_vertices[0] = a;
	initial_vertices[1] = b;
	initial_vertices[2] = c;
	initial_vertices[3] = d;
	return 0;
}

static int int_compare(const void* va, const void* vb)
{
	return *(int*)va - *(int*)vb;
}

static int initialize_convex_hull(int num_points, const double (*points)[4], int8_t facets[][3], double plane_normal[][4], double* barycentre, int* initial_vertices, bool* processed, int* p_num_facets)
{
	int ret = initial_simplex(num_points, points, initial_vertices);
	assert(ret == 0);

	int a = initial_vertices[0];
	int b = initial_vertices[1];
	int c = initial_vertices[2];
	int d = initial_vertices[3];
	int e = -1;

	if (!find_fifth_point(num_points, points, a, b, c, d, &e))
		return -1;

	double _centre[3];
	//assert(circumsphere_centre((double*)points[b], (double*)points[c], (double*)points[e], _centre));
	//assert(circumsphere_centre((double*)points[b], (double*)points[d], (double*)points[e], _centre));
	//assert(circumsphere_centre((double*)points[c], (double*)points[d], (double*)points[e], _centre));

	if (!circumsphere_centre((double*)points[b], (double*)points[c], (double*)points[e], _centre)) return -101;
	if (!circumsphere_centre((double*)points[b], (double*)points[d], (double*)points[e], _centre)) return -101;
	if (!circumsphere_centre((double*)points[c], (double*)points[d], (double*)points[e], _centre)) return -101;

	initial_vertices[4] = e;

#ifdef DEBUG
	printf("initial vertices: %d %d %d %d %d\n", a, b, c, d, e);
#endif


	qsort(initial_vertices, 5, sizeof(int), int_compare);

	b = initial_vertices[1];
	c = initial_vertices[2];
	d = initial_vertices[3];
	e = initial_vertices[4];

	for (int i=0;i<5;i++)
	{
		processed[initial_vertices[i]] = true;
		barycentre[0] += points[initial_vertices[i]][0];
		barycentre[1] += points[initial_vertices[i]][1];
		barycentre[2] += points[initial_vertices[i]][2];
		barycentre[3] += points[initial_vertices[i]][3];
	}
	barycentre[0] /= 5;
	barycentre[1] /= 5;
	barycentre[2] /= 5;
	barycentre[3] /= 5;

	int num_facets = 0;
	_add_facet(points, b, c, d, facets[num_facets], plane_normal[num_facets], barycentre); num_facets++;
	_add_facet(points, b, c, e, facets[num_facets], plane_normal[num_facets], barycentre); num_facets++;
	_add_facet(points, b, d, e, facets[num_facets], plane_normal[num_facets], barycentre); num_facets++;
	_add_facet(points, c, d, e, facets[num_facets], plane_normal[num_facets], barycentre); num_facets++;

	*p_num_facets = num_facets;
	return 0;
}

static int get_convex_hull_4D(int num_points, const double (*points)[4], int* p_num_facets, int8_t facets[][3])
{
	assert(num_points <= MAX_POINTS);
	double plane_normal[MAXF][4];
	double barycentre[4] = {0, 0, 0, 0};

	int num_facets = 0;
	int initial_vertices[5];
	bool processed[MAX_POINTS] = {false};
	int ret = initialize_convex_hull(num_points, (const double (*)[4])points, facets, plane_normal, barycentre, initial_vertices, processed, &num_facets);
	if (ret != 0)
		return ret;
	//assert(ret == 0);
	assert(num_facets >= 2);

	int8_t to_add[MAXF][3];
	uint8_t edgehit[MAX_POINTS][MAX_POINTS];
#define VISIBLE 1
#define INVISIBLE 2
#define BOTH 3

	for (int i = 0;i<num_points;i++)
	{
		if (processed[i])
			continue;
		processed[i] = true;

		int num_to_add = 0;
		memset(edgehit, 0, MAX_POINTS*MAX_POINTS * sizeof(uint8_t));

		for (int j = 0;j<num_facets;j++)
		{
			int b = facets[j][0];
			int c = facets[j][1];
			int d = facets[j][2];

			int vertices[3][2] = {	{ b, c},
						{ b, d},
						{ c, d}	};

			bool both[3] = {false, false, false};

			double distance = point_plane_distance_4D(points[i], plane_normal[j]);
			bool vis = distance > TOLERANCE;
#ifdef DEBUG
printf("point %d facet %d %d %d %d distance %f\n", i, 0, b, c, d, distance);
#endif
			if (vis)
			{
				for (int k=0;k<3;k++)
				{
					int u = vertices[k][0];
					int v = vertices[k][1];
					edgehit[u][v] |= VISIBLE;
					both[k] = edgehit[u][v] == BOTH;
				}

#ifdef DEBUG
				printf("removing: %d %d %d %d\n", 0, b, c, d);
#endif

				memcpy(facets[j], facets[num_facets-1], 3 * sizeof(int8_t));
				memcpy(plane_normal[j], plane_normal[num_facets-1], 4 * sizeof(double));
				num_facets--;
				j--;
			}
			else
			{
				for (int k=0;k<3;k++)
				{
					int u = vertices[k][0];
					int v = vertices[k][1];
					edgehit[u][v] |= INVISIBLE;
					both[k] = edgehit[u][v] == BOTH;
				}
			}

			for (int k=0;k<3;k++)
			{
				if (both[k])
				{
					if (i > vertices[k][1])
					{
						to_add[num_to_add][0] = vertices[k][0];
						to_add[num_to_add][1] = vertices[k][1];
						to_add[num_to_add][2] = i;
					}
					else if (i > vertices[k][0])
					{
						to_add[num_to_add][0] = vertices[k][0];
						to_add[num_to_add][1] = i;
						to_add[num_to_add][2] = vertices[k][1];
					}
					else
					{
						to_add[num_to_add][0] = i;
						to_add[num_to_add][1] = vertices[k][0];
						to_add[num_to_add][2] = vertices[k][1];

					}

#ifdef DEBUG
					double _centre[3];
					if (!circumsphere_centre((double*)points[to_add[num_to_add][0]], (double*)points[to_add[num_to_add][1]], (double*)points[to_add[num_to_add][2]], _centre))
					{
							printf("coplanar facet added! %d %d %d %d\n", 0, to_add[num_to_add][0], to_add[num_to_add][1], to_add[num_to_add][2]);
							return -1;
					}
#endif

					num_to_add++;
					//assert(num_to_add < MAXF);
					if (num_to_add >= MAXF)
					{
#ifdef DEBUG
						printf("number of facets to add exceeds maximum\n");
#endif
						return -1;
					}
				}
			}
		}

		for (int j = 0;j<num_to_add;j++)
		{
			//assert(num_facets < MAXF);
			if (num_facets >= MAXF)
			{
#ifdef DEBUG
				printf("number of facets added exceeds maximum\n");
#endif
				return -1;
			}

			_add_facet((const double (*)[4])points, to_add[j][0], to_add[j][1], to_add[j][2], facets[num_facets], plane_normal[num_facets], barycentre);
			num_facets++;
		}
	}


	*p_num_facets = num_facets;
	if (num_facets == 0)
	{
#ifdef DEBUG
		printf("zero facets!\n");
#endif
		return -1;
	}
	return ret;
}

static double triangle_area(double* x1, double* x2, double* x3)
{
	//http://mathworld.wolfram.com/TriangleArea.html

	double u[3] = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
	double v[3] = {x3[0] - x2[0], x3[1] - x2[1], x3[2] - x2[2]};

	double c[3];
	cross_product(u, v, c);
	return sqrt(norm_squared(c)) / 2;
}

static int calculate_voronoi_areas(int num_points, const double (*points)[4], int num_facets, int8_t facets[][3], double* areas)
{
#ifdef DEBUG
	printf("num facets: %d\n", num_facets);
#endif

	int ret = 0;
	bool point_hit[MAX_POINTS] = {false};
	double circumcentres[MAXF][3];
	for (int i=0;i<num_facets;i++)
	{
		int b = facets[i][0];
		int c = facets[i][1];
		int d = facets[i][2];

		point_hit[b] = true;
		point_hit[c] = true;
		point_hit[d] = true;

#ifdef DEBUG
		printf("facet %d: %d %d %d %d\n", i, 0, b, c, d);
#endif

		if (!circumsphere_centre((double*)points[b], (double*)points[c], (double*)points[d], circumcentres[i]))
		{
#ifdef DEBUG
			printf("planar facet: %d %d %d %d\n", 0, b, c, d);
#endif
			return -1;
		}
	}

	int8_t adj[MAX_POINTS][MAX_POINTS];
	memset(adj, -1, sizeof(int8_t) * MAX_POINTS * MAX_POINTS);

	int8_t adjf[MAX_POINTS][MAX_POINTS];
	memset(adjf, -1, sizeof(int8_t) * MAX_POINTS * MAX_POINTS);

	int8_t first_facet[MAX_POINTS];
	memset(first_facet, -1, sizeof(int8_t) * MAX_POINTS);

	for (int i=0;i<num_facets;i++)
	{
		int b = facets[i][0];
		int c = facets[i][1];
		int d = facets[i][2];

		if (first_facet[b] == -1)	first_facet[b] = i;
		if (first_facet[c] == -1)	first_facet[c] = i;
		if (first_facet[d] == -1)	first_facet[d] = i;

		if (adjf[b][c] == -1)
		{
			adjf[b][c] = i;
			adj[b][c] = d;
		}
		else if (adjf[c][b] == -1)
		{
			adjf[c][b] = i;
			adj[c][b] = d;
		}
		else
			return -1;

		if (adjf[c][d] == -1)
		{
			adjf[c][d] = i;
			adj[c][d] = b;
		}
		else if (adjf[d][c] == -1)
		{
			adjf[d][c] = i;
			adj[d][c] = b;
		}
		else
			return -1;

		if (adjf[d][b] == -1)
		{
			adjf[d][b] = i;
			adj[d][b] = c;
		}
		else if (adjf[b][d] == -1)
		{
			adjf[b][d] = i;
			adj[b][d] = c;
		}
		else
			return -1;
	}

	areas[0] = INFINITY;
	for (int i=1;i<num_points;i++)
	{
		if (!point_hit[i])
		{
			areas[i] = 0;
			continue;
		}

		int8_t localfacets[MAXF];
		int pivot = facets[first_facet[i]][0];
		if (pivot == i)
			pivot = facets[first_facet[i]][1];

		int prev = facets[first_facet[i]][0];
		if (prev == i || prev == pivot)
			prev = facets[first_facet[i]][1];
		if (prev == i || prev == pivot)
			prev = facets[first_facet[i]][2];

		int orig = prev;
		int num_local = 1;

		localfacets[0] = first_facet[i];
		for (int j=1;pivot != orig;j++, num_local++)
		{
			if (adj[i][pivot] != -1 && adj[i][pivot] != prev)
			{
				localfacets[j] = adjf[i][pivot];
				prev = pivot;
				pivot = adj[i][pivot];
			}
			else if (adj[pivot][i] != -1 && adj[pivot][i] != prev)
			{
				localfacets[j] = adjf[pivot][i];
				prev = pivot;
				pivot = adj[pivot][i];
			}
			else
				return -1;
		}

		double area = 0.0;
		int u = localfacets[0];
		int v = localfacets[1];
		for (int j=2;j<num_local;j++)
		{
			int w = localfacets[j];
			area += triangle_area(circumcentres[u], circumcentres[v], circumcentres[w]);
			v = w;
		}

		areas[i] = area;
	}

	return ret;
}

typedef struct
{
	double area;
	double dist;
	int index;
} sorthelper_t;

static int sorthelper_compare(const void* va, const void* vb)
{
	sorthelper_t* a = (sorthelper_t*)va;
	sorthelper_t* b = (sorthelper_t*)vb;

	if (a->area > b->area)
		return -1;

	if (a->area < b->area)
		return 1;

	if (a->dist < b->dist)
		return -1;

	return 1;
}

int calculate_neighbour_ordering(int num_points, const double (*_points)[3], int8_t* ordering)
{
assert(num_points <= MAX_POINTS);

	double points[MAX_POINTS][4];
	for (int i = 0;i<num_points;i++)
	{
		double x = _points[i][0] - _points[0][0];
		double y = _points[i][1] - _points[0][1];
		double z = _points[i][2] - _points[0][2];
		points[i][0] = x;
		points[i][1] = y;
		points[i][2] = z;
		points[i][3] = x*x + y*y + z*z;//norm_squared((double*)_points[i]);

#ifdef DEBUG
		printf("point %d: %f\t%f\t%f\t%f\n", i, x, y, z, x*x + y*y + z*z);
#endif
	}

	int num_facets = 0;
	int8_t simplex[MAXF][3];
	int ret = get_convex_hull_4D(num_points, (const double (*)[4])points, &num_facets, simplex);
	//assert(ret == 0);
	if (ret != 0)
		return ret;	//todo: replace with assert again

	double areas[MAX_POINTS];
	ret = calculate_voronoi_areas(num_points, (const double (*)[4])points, num_facets, simplex, areas);
	if (ret != 0)
		return ret;

	sorthelper_t data[MAX_POINTS];
	for (int i=0;i<num_points;i++)
	{
		assert(areas[i] == areas[i]);
		//printf("%d: %f\n", i, areas[i]);

		data[i].area = areas[i];
		data[i].dist = norm_squared((double*)_points[i]);
		data[i].index = i;
	}

	qsort(data, num_points, sizeof(sorthelper_t), sorthelper_compare);

#ifdef DEBUG
	for (int i=0;i<num_points;i++)
		printf("%d %f\n", data[i].index, data[i].area);
#endif

	for (int i=0;i<num_points;i++)
		ordering[i] = data[i].index;

	return ret;
}

