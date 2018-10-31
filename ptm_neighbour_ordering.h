#ifndef PTM_NEIGHBOUR_ORDERING_H
#define PTM_NEIGHBOUR_ORDERING_H

namespace ptm {

int calculate_neighbour_ordering(void* voronoi_handle, int num_points, const double (*_points)[3], int8_t* ordering);

int calculate_diamond_neighbour_ordering(	void* _voronoi_handle, size_t atom_index, int (get_neighbours)(void* vdata, size_t atom_index, int num, size_t* nbr_indices, int32_t* numbers, double (*nbr_pos)[3]), void* nbrlist,
						int8_t* ordering, double (*points)[3], int32_t* numbers);

void* voronoi_initialize_local();
void voronoi_uninitialize_local(void* ptr);

}

#endif

