#ifndef NEIGHBOUR_ORDERING_HPP
#define NEIGHBOUR_ORDERING_HPP

int calculate_neighbour_ordering(void* voronoi_handle, int num_points, const double (*_points)[3], int8_t* ordering);

void* voronoi_initialize_local();
void voronoi_uninitialize_local(void* ptr);

#endif

