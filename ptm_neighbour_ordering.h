#ifndef PTM_NEIGHBOUR_ORDERING_H
#define PTM_NEIGHBOUR_ORDERING_H

#include <cstddef>

namespace ptm {

int calculate_neighbour_ordering(	void* _voronoi_handle, size_t atom_index, int min_points, int (get_neighbours)(void* vdata, size_t atom_index, int num, size_t* nbr_indices, int32_t* numbers,
					double (*nbr_pos)[3]), void* nbrlist, bool calc_solid_angles,
					size_t* nbr_indices, double (*points)[3], int32_t* numbers);

int calculate_two_shell_neighbour_ordering(	void* _voronoi_handle, size_t atom_index, int (get_neighbours)(void* vdata, size_t atom_index, int num, size_t* nbr_indices, int32_t* numbers, double (*nbr_pos)[3]), void* nbrlist,
						int num_inner, int num_outer, int max_snbrs, bool calc_solid_angles,
						size_t* nbr_indices, double (*points)[3], int32_t* numbers);

void* voronoi_initialize_local();
void voronoi_uninitialize_local(void* ptr);

}

#endif

