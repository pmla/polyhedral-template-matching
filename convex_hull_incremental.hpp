#ifndef CONVEX_HULL_INCREMENTAL_HPP
#define CONVEX_HULL_INCREMENTAL_HPP


#include <cstdint>
#include <cstdbool>

#define MAXF 24
typedef struct
{
	int8_t facets[MAXF][3];
	double plane_normal[MAXF][3];
	bool processed[15];
	int initial_vertices[4];
	double barycentre[3];
	int num_facets;
	int num_prev;
	bool ok;

} convexhull_t;

void add_facet(const double (*points)[3], int a, int b, int c, int8_t* facet, double* plane_normal, double* barycentre);
int get_convex_hull(int num_points, const double (*points)[3], int num_expected_facets, convexhull_t* ch, int8_t simplex[][3]);

#endif

