#ifndef CONVEX_HULL_HPP
#define CONVEX_HULL_HPP

void add_facet(const double* points, int a, int b, int c, int8_t* facet, double* plane_normal, double* barycentre);
bool get_convex_hull(int num_points, double* points, int num_expected_facets, int8_t facet[][3]);

#endif

