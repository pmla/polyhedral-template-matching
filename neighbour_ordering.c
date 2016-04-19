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


/* int triangular(int n)
{
	return n*(n+1)/2;
}

static int tetrahedral(int n)
{
	return n*(n+1)*(n+2) / 6;
}

static int triple_index(int i, int j, int k)
{
	return tetrahedral(MAX_POINTS - 2) - tetrahedral(MAX_POINTS - i - 2) + (triangular(MAX_POINTS-i-2) - triangular(MAX_POINTS-1-j)) + k-j-1;
}

static void setb(uint8_t* buf, int i)
{
	int wordpos = i / 8;
	int bitpos = i % 8;

	buf[wordpos] |= 1 << bitpos;
}

static int getb(uint8_t* buf, int i)
{
	int wordpos = i / 8;
	int bitpos = i % 8;
	return (buf[wordpos] >> bitpos) & 1;
}

static double vector_distance_2(double* a, double* b)
{
	double x = a[0] - b[0];
	double y = a[1] - b[1];
	double z = a[2] - b[2];

	return x*x + y*y + z*z;
}
*/

static uint16_t index_i[MAX_POINTS] = {153,289,409,514,605,683,749,804,849,885,913,934,949,959,965,968,969,969,969};
static uint8_t index_j[MAX_POINTS] = {172,155,139,124,110,97,85,74,64,55,47,40,34,29,25,22,20,19,19};

static int triple_index(int i, int j, int k)
{
	return index_i[i] - index_j[j] + k;
}

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

double dethelper(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4)
{
	return x1*y2*z3-x1*y2*z4-x1*y3*z2+x1*y3*z4+x1*y4*z2-x1*y4*z3
		-x2*y1*z3+x2*y1*z4+x2*y3*z1-x2*y3*z4-x2*y4*z1+x2*y4*z3
		+x3*y1*z2-x3*y1*z4-x3*y2*z1+x3*y2*z4+x3*y4*z1-x3*y4*z2
		-x4*y1*z2+x4*y1*z3+x4*y2*z1-x4*y2*z3-x4*y3*z1+x4*y3*z2;
}

static void circumsphere_centre(double* a, double* b, double* c, double* d, double* centre)
{
	//formula adapted from http://mathworld.wolfram.com/Circumsphere.html

	double det = dethelper(a[0], b[0], c[0], d[0], a[1], b[1], c[1], d[1], a[2], b[2], c[2], d[2]);

	//double sa = norm_squared(a);
	//double sb = norm_squared(b);
	//double sc = norm_squared(c);
	//double sd = norm_squared(d);

	//norm squared is stored in fourth element
	double sa = a[3];
	double sb = b[3];
	double sc = c[3];
	double sd = d[3];

	double detx = dethelper(sa, sb, sc, sd, a[1], b[1], c[1], d[1], a[2], b[2], c[2], d[2]);
	double dety = dethelper(sa, sb, sc, sd, a[0], b[0], c[0], d[0], a[2], b[2], c[2], d[2]);
	double detz = dethelper(sa, sb, sc, sd, a[0], b[0], c[0], d[0], a[1], b[1], c[1], d[1]);

	centre[0] = detx / (2 * det);
	centre[1] = -dety / (2 * det);
	centre[2] = detz / (2 * det);
}

static void cross_product_4D(double* r0, double* r1, double* r2, double* p)
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

static void calculate_plane_normal_4D(const double (*points)[4], int a, int b, int c, int d, double* plane_normal)
{
	double u[4] = {	points[b][0] - points[a][0],
			points[b][1] - points[a][1],
			points[b][2] - points[a][2],
			points[b][3] - points[a][3]	};

	double v[4] = {	points[c][0] - points[a][0],
			points[c][1] - points[a][1],
			points[c][2] - points[a][2],
			points[c][3] - points[a][3]	};

	double w[4] = {	points[d][0] - points[a][0],
			points[d][1] - points[a][1],
			points[d][2] - points[a][2],
			points[d][3] - points[a][3]	};

	cross_product_4D(u, v, w, plane_normal);
	double norm = sqrt(norm_squared_4D(plane_normal));
	plane_normal[0] /= norm;
	plane_normal[1] /= norm;
	plane_normal[2] /= norm;
	plane_normal[3] /= norm;
}

static double point_plane_distance_4D(const double* w, const double* plane_point, const double* plane_cross)
{
	return	  plane_cross[0] * (plane_point[0] - w[0])
		+ plane_cross[1] * (plane_point[1] - w[1])
		+ plane_cross[2] * (plane_point[2] - w[2])
		+ plane_cross[3] * (plane_point[3] - w[3]);
}

static bool visible(const double* w, const double* plane_point, const double* plane_normal)
{
	return point_plane_distance_4D(w, plane_point, plane_normal) > TOLERANCE;
}

static void _add_facet(const double (*points)[4], int a, int b, int c, int d, int8_t* facet, double* plane_normal, double* barycentre)
{
#ifdef DEBUG
	printf("adding: %d %d %d %d\n", a, b, c, d);
	assert(a < b);
	assert(a < c);
	assert(a < d);
	assert(b < c);
	assert(b < d);
	assert(c < d);
#endif

	facet[0] = a;
	facet[1] = b;
	facet[2] = c;
	facet[3] = d;

	calculate_plane_normal_4D(points, a, b, c, d, plane_normal);
	if (visible(barycentre, points[a], plane_normal))
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

static int initialize_convex_hull(int num_points, const double (*points)[4], int8_t facets[][4], double plane_normal[][4], double* barycentre, int* initial_vertices, bool* processed, int* p_num_facets)
{
	int ret = initial_simplex(num_points, points, initial_vertices);
	assert(ret == 0);

	int a = initial_vertices[0];
	int b = initial_vertices[1];
	int c = initial_vertices[2];
	int d = initial_vertices[3];
	int e = -1;

	for (int i=0;i<num_points;i++)
	{
		if (i != a &&i != b &&i != c &&i != d)
		{
			e = i;
			break;
		}
	}
	initial_vertices[4] = e;
	qsort(initial_vertices, 5, sizeof(int), int_compare);

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
	for (int i=0;i<5;i++)
	{
		int vertices[4];
		int n=0;
		for (int j=0;j<5;j++)
			if (j != i)
				vertices[n++] = initial_vertices[j];

		_add_facet(points, vertices[0], vertices[1], vertices[2], vertices[3], facets[num_facets], plane_normal[num_facets], barycentre);
		num_facets++;
	}

	*p_num_facets = num_facets;
	return 0;
}

static int get_convex_hull_4D(int num_points, const double (*points)[4], int* p_num_facets, int8_t facets[][4])
{
	assert(num_points <= MAX_POINTS);
	double plane_normal[MAXF][4];
	double barycentre[4] = {0, 0, 0, 0};

	int num_facets = 0;
	int initial_vertices[5];
	bool processed[MAX_POINTS] = {false};
	int ret = initialize_convex_hull(num_points, (const double (*)[4])points, facets, plane_normal, barycentre, initial_vertices, processed, &num_facets);
	assert(ret == 0);
	assert(num_facets >= 2);

	int8_t to_add[MAXF][4];
	//uint8_t edge_visible[122];
	//uint8_t edge_invisible[122];
	uint8_t edgehit[1000];
#define VISIBLE 1
#define INVISIBLE 2
#define BOTH 3

	for (int i = 0;i<num_points;i++)
	{
		if (processed[i])
			continue;
		processed[i] = true;

		int num_to_add = 0;
		//memset(edge_visible, 0, sizeof(uint8_t) * 122);
		//memset(edge_invisible, 0, sizeof(uint8_t) * 122);
		memset(edgehit, 0, 1000 * sizeof(uint8_t));

		for (int j = 0;j<num_facets;j++)
		{
			int a = facets[j][0];
			int b = facets[j][1];
			int c = facets[j][2];
			int d = facets[j][3];

			int vertices[4][3] = {	{ a, b, c},
						{ a, b, d},
						{ a, c, d},
						{ b, c, d}	};

			bool both[4] = {false, false, false, false};

			int ts[4];
			for (int k=0;k<4;k++)
				ts[k] = triple_index(vertices[k][0], vertices[k][1], vertices[k][2]);

			double distance = point_plane_distance_4D(points[i], points[a], plane_normal[j]);
			bool vis = distance > TOLERANCE;
			if (vis)
			{
				for (int k=0;k<4;k++)
				{
					//setb(edge_visible, ts[k]);
					//both[k] = getb(edge_invisible, ts[k]);

					edgehit[ts[k]] |= VISIBLE;
					both[k] = edgehit[ts[k]] == BOTH;
				}

#ifdef DEBUG
				printf("removing: %d %d %d %d\n", a, b, c, d);
#endif

				memcpy(facets[j], facets[num_facets-1], 4 * sizeof(int8_t));
				memcpy(plane_normal[j], plane_normal[num_facets-1], 4 * sizeof(double));
				num_facets--;
				j--;
			}
			else
			{
				for (int k=0;k<4;k++)
				{
					//setb(edge_invisible, ts[k]);
					//both[k] = getb(edge_visible, ts[k]);

					edgehit[ts[k]] |= INVISIBLE;
					both[k] = edgehit[ts[k]] == BOTH;
				}
			}

			for (int k=0;k<4;k++)
			{
				if (both[k] && vertices[k][0] == 0)
				{
					if (i > vertices[k][2])
					{
						to_add[num_to_add][0] = vertices[k][0];
						to_add[num_to_add][1] = vertices[k][1];
						to_add[num_to_add][2] = vertices[k][2];
						to_add[num_to_add][3] = i;
					}
					else if (i > vertices[k][1])
					{
						to_add[num_to_add][0] = vertices[k][0];
						to_add[num_to_add][1] = vertices[k][1];
						to_add[num_to_add][2] = i;
						to_add[num_to_add][3] = vertices[k][2];
					}
					else
					{
						to_add[num_to_add][0] = vertices[k][0];
						to_add[num_to_add][1] = i;
						to_add[num_to_add][2] = vertices[k][1];
						to_add[num_to_add][3] = vertices[k][2];
					}

					num_to_add++;
					assert(num_to_add < MAXF);
				}
			}
		}

		for (int j = 0;j<num_to_add;j++)
		{
			assert(num_facets < MAXF);

			//if (to_add[j][0] == 0)
			{
				_add_facet((const double (*)[4])points, to_add[j][0], to_add[j][1], to_add[j][2], to_add[j][3], facets[num_facets], plane_normal[num_facets], barycentre);
				num_facets++;
			}
		}
	}

	for (int i=0;i<num_facets;i++)
	{
		if (facets[i][0] != 0)
		{
			memcpy(&facets[i], &facets[num_facets - 1], 4 * sizeof(int8_t));
			num_facets--;
			i--;
		}
	}

	*p_num_facets = num_facets;
	return ret;
}

static double triangle_area(double* x1, double* x2, double* x3)
{
	//http://mathworld.wolfram.com/TriangleArea.html

	double u[3] = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
	double v[3] = {x3[0] - x2[0], x3[1] - x2[1], x3[2] - x2[2]};

	double c[3];
	cross_product(u, v, c);
	return sqrt(fabs(norm_squared(c))) / 2;	//fabs guards against very small negative numbers
}

static int calculate_voronoi_areas(int num_points, const double (*points)[4], int num_facets, int8_t facets[][4], double* areas)
{
	int ret = 0;
	bool point_hit[MAX_POINTS] = {false};
	double circumcentres[MAXF][3];
	for (int i=0;i<num_facets;i++)
	{
		int a = facets[i][0];
		int b = facets[i][1];
		int c = facets[i][2];
		int d = facets[i][3];

		point_hit[a] = true;
		point_hit[b] = true;
		point_hit[c] = true;
		point_hit[d] = true;

		circumsphere_centre((double*)points[a], (double*)points[b], (double*)points[c], (double*)points[d], circumcentres[i]);
	}

	areas[0] = INFINITY;
	for (int i=1;i<num_points;i++)
	{
		if (!point_hit[i])
		{
			areas[i] = 0;
			continue;
		}

		int num_local = 0;
		int8_t localfacets[MAXF];
		for (int j=0;j<num_facets;j++)
		{
			int a = facets[j][0];
			int b = facets[j][1];
			int c = facets[j][2];
			int d = facets[j][3];

			assert(a == 0);
			if (b == i || c == i || d == i)
				localfacets[num_local++] = j;
		}

		int _b = i;
		int _c = facets[localfacets[0]][1];
		if (_c == _b)
			_c = facets[localfacets[0]][2];
		assert(_c != _b);

		for (int j=1;j<num_local;j++)
		{
			bool found = false;
			for (int k=j;k<num_local;k++)
			{

				int a = facets[localfacets[k]][0];
				int b = facets[localfacets[k]][1];
				int c = facets[localfacets[k]][2];
				int d = facets[localfacets[k]][3];

				assert(a == 0);
				if ((b == _b && c == _c) || (b == _c && c == _b))
				{
					int8_t temp = localfacets[j];
					localfacets[j] = localfacets[k];
					localfacets[k] = temp;

					found = true;
					_c = d;
					break;
				}
				else if ((b == _b && d == _c) || (b == _c && d == _b))
				{
					int8_t temp = localfacets[j];
					localfacets[j] = localfacets[k];
					localfacets[k] = temp;

					found = true;
					_c = c;
					break;
				}
				else if ((c == _b && d == _c) || (c == _c && d == _b))
				{
					int8_t temp = localfacets[j];
					localfacets[j] = localfacets[k];
					localfacets[k] = temp;

					found = true;
					_c = b;
					break;
				}
			}

			if (!found) return -1;
			//this assertion will probably fail if atom lies at the surface (unbounded voronoi cell)
			assert(found);
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
		double x = _points[i][0];
		double y = _points[i][1];
		double z = _points[i][2];
		points[i][0] = x;
		points[i][1] = y;
		points[i][2] = z;
		points[i][3] = norm_squared((double*)_points[i]);

#ifdef DEBUG
		printf("point %d: %f\t%f\t%f\t%f\n", i, x, y, z, x*x + y*y + z*z);
#endif
	}

	int num_facets = 0;
	int8_t simplex[MAXF][4];
	int ret = get_convex_hull_4D(num_points, (const double (*)[4])points, &num_facets, simplex);
	//assert(ret == 0);
	if (ret != 0)
		return -1;	//todo: replace with assert again

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

