#include <cstdint>
#include <cassert>
#include <cstring>
#include <cfloat>
#include <cmath>

#define VISIBLE 1
#define INVISIBLE 2
#define BOTH 3
#define MAXF 24


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

static void calculate_plane_normal(const double* points, int a, int b, int c, double* plane_normal)
{
	double u[3] = {	points[b * 3 + 0] - points[a * 3 + 0],
			points[b * 3 + 1] - points[a * 3 + 1],
			points[b * 3 + 2] - points[a * 3 + 2]	};

	double v[3] = {	points[c * 3 + 0] - points[a * 3 + 0],
			points[c * 3 + 1] - points[a * 3 + 1],
			points[c * 3 + 2] - points[a * 3 + 2]	};

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

static bool calc_max_extent(int num_points, double* points, int* min_index, int* max_index)
{
	for (int j=0;j<3;j++)
	{
		double dmin = DBL_MAX, dmax = -DBL_MAX;
		int imin = 0, imax = 0;

		for (int i = 0;i<num_points;i++)
		{
			double d = points[i * 3 + j];
			if (d < dmin)
			{
				dmin = d;
				imin = i;
			}
			if (d > dmax)
			{
				dmax = d;
				imax = i;
			}
		}

		if (imin == imax)
			return false;	//degenerate point set

		min_index[j] = imin;
		max_index[j] = imax;
	}

	return true;
}

static bool find_third_point(int num_points, double* points, int a, int b, double tolerance, int* p_c)
{
	double* x1 = &points[a * 3];
	double* x2 = &points[b * 3];

	double x2x1[3] = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
	double ns_x2x1 = norm_squared(x2x1);

	int bi = -1;
	double max_dist = 0.0;
	for (int i = 0;i<num_points;i++)
	{
		if (i == a || i == b)
			continue;

		double* x0 = &points[i * 3];

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
	return max_dist > tolerance;
}

static bool find_fourth_point(int num_points, double* points, int a, int b, int c, double tolerance, int* p_d)
{
	double plane_normal[3];
	calculate_plane_normal(points, a, b, c, plane_normal);


	int bi = -1;
	double max_dist = 0.0;
	for (int i = 0;i<num_points;i++)
	{
		if (i == a || i == b || i == c)
			continue;

		double* x0 = &points[i * 3];
		double dist = fabs(point_plane_distance(x0, &points[a*3 + 0], plane_normal));
		if (dist > max_dist)
		{
			max_dist = dist;
			bi = i;
		}
	}

	*p_d = bi;
	return max_dist > tolerance;
}

static bool initial_simplex(int num_points, double* points, double tolerance, int* initial_vertices)
{
	int min_index[3] = {0};
	int max_index[3] = {0};
	if (!calc_max_extent(num_points, points, min_index, max_index))
		return false;

	int bi = -1;
	double max_dist = 0.0;
	for (int i = 0;i<3;i++)
	{
		int a = min_index[i], b = max_index[i];
		double delta[3] = {	points[a*3 + 0] - points[b*3 + 0],
					points[a*3 + 1] - points[b*3 + 1],
					points[a*3 + 2] - points[b*3 + 2]	};
		double dist = norm_squared(delta);
		if (dist > max_dist)
		{
			bi = i;
			max_dist = dist;
		}
	}

	//first two points are (a, b)
	int a = min_index[bi], b = max_index[bi], c = -1, d = -1;

	if (!find_third_point(num_points, points, a, b, tolerance, &c))
		return false;

	if (!find_fourth_point(num_points, points, a, b, c, tolerance, &d))
		return false;

	initial_vertices[0] = a;
	initial_vertices[1] = b;
	initial_vertices[2] = c;
	initial_vertices[3] = d;
	return true;
}

static bool visible(const double* w, const double* plane_point, const double* plane_normal)
{
	return point_plane_distance(w, plane_point, plane_normal) > 0;
}

void add_facet(const double* points, int a, int b, int c, int8_t* facet, double* plane_normal, double* barycentre)
{
	calculate_plane_normal(points, a, b, c, plane_normal);
	if (visible(barycentre, &points[a * 3], plane_normal))
	{
		plane_normal[0] = -plane_normal[0];
		plane_normal[1] = -plane_normal[1];
		plane_normal[2] = -plane_normal[2];

		facet[0] = b;
		facet[1] = a;
		facet[2] = c;
	}
	else
	{
		facet[0] = a;
		facet[1] = b;
		facet[2] = c;
	}
}

int get_convex_hull(int num_points, double* points, int num_expected_facets, int8_t simplex[][3])
{
	assert(num_points == 7 || num_points == 13 || num_points == 15);

	double tolerance = 1E-8;
	int8_t facets[MAXF][3];
	double plane_normal[MAXF][3];
	int8_t edge_visible[15][15];
	bool processed[15] = {false};

	int num_to_add = 0;
	int8_t to_add[MAXF][3];

	double barycentre[3] = {0};
	int initial_vertices[4] = {0, 1, 2, num_points - 1};
	if (!initial_simplex(num_points, points, tolerance, initial_vertices))
		return false;

	for (int i = 0;i<4;i++)
	{
		int a = initial_vertices[i];
		processed[a] = true;

		barycentre[0] += points[a * 3 + 0];
		barycentre[1] += points[a * 3 + 1];
		barycentre[2] += points[a * 3 + 2];
	}
	barycentre[0] /= 4;
	barycentre[1] /= 4;
	barycentre[2] /= 4;

	int num_facets = 0;
	add_facet(points, initial_vertices[0], initial_vertices[1], initial_vertices[2], facets[num_facets], plane_normal[num_facets], barycentre); num_facets++;
	add_facet(points, initial_vertices[0], initial_vertices[1], initial_vertices[3], facets[num_facets], plane_normal[num_facets], barycentre); num_facets++;
	add_facet(points, initial_vertices[0], initial_vertices[2], initial_vertices[3], facets[num_facets], plane_normal[num_facets], barycentre); num_facets++;
	add_facet(points, initial_vertices[1], initial_vertices[2], initial_vertices[3], facets[num_facets], plane_normal[num_facets], barycentre); num_facets++;

	for (int i = 0;i<num_points;i++)
	{
		if (processed[i])
			continue;
		processed[i] = true;

		num_to_add = 0;
		memset(edge_visible, 0, sizeof(int8_t) * 15 * 15);
		for (int j = 0;j<num_facets;j++)
		{
			int a = facets[j][0];
			int b = facets[j][1];
			int c = facets[j][2];

			int u = 0, v = 0, w = 0;

			double distance = point_plane_distance(&points[i * 3], &points[a * 3], plane_normal[j]);
			bool vis = distance > tolerance;
			if (vis)
			{
				u = edge_visible[a][b] |= VISIBLE;
				edge_visible[b][a] |= VISIBLE;

				v = edge_visible[b][c] |= VISIBLE;
				edge_visible[c][b] |= VISIBLE;

				w = edge_visible[c][a] |= VISIBLE;
				edge_visible[a][c] |= VISIBLE;

				memcpy(facets[j], facets[num_facets-1], 3 * sizeof(int8_t));
				memcpy(plane_normal[j], plane_normal[num_facets-1], 3 * sizeof(double));
				num_facets--;
				j--;
			}
			else
			{
				u = edge_visible[a][b] |= INVISIBLE;
				edge_visible[b][a] |= INVISIBLE;

				v = edge_visible[b][c] |= INVISIBLE;
				edge_visible[c][b] |= INVISIBLE;

				w = edge_visible[c][a] |= INVISIBLE;
				edge_visible[a][c] |= INVISIBLE;
			}

			if (u == BOTH)
			{
				to_add[num_to_add][0] = i;
				to_add[num_to_add][1] = a;
				to_add[num_to_add][2] = b;
				num_to_add++;
			}

			if (v == BOTH)
			{
				to_add[num_to_add][0] = i;
				to_add[num_to_add][1] = b;
				to_add[num_to_add][2] = c;
				num_to_add++;
			}

			if (w == BOTH)
			{
				to_add[num_to_add][0] = i;
				to_add[num_to_add][1] = c;
				to_add[num_to_add][2] = a;
				num_to_add++;
			}
		}

		for (int j = 0;j<num_to_add;j++)
		{

			if (num_facets >= MAXF)
				return false;

			add_facet(points, to_add[j][0], to_add[j][1], to_add[j][2], facets[num_facets], plane_normal[num_facets], barycentre); num_facets++;
		}
	}


	if (num_facets != num_expected_facets)
		return false;			//incorrect number of facets in convex hull

	for (int i=0;i<num_facets;i++)
	{
		int a = facets[i][0];
		int b = facets[i][1];
		int c = facets[i][2];
		if (a == num_points - 1 || b == num_points - 1 || c == num_points - 1)
			return false;		//central atom contained in convex hull

		simplex[i][0] = a;
		simplex[i][1] = b;
		simplex[i][2] = c;
	}

	return true;
}

