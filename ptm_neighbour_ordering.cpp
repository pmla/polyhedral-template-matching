//todo: normalize vertices

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <unordered_set>
#include "ptm_constants.h"
#include "ptm_voronoi_cell.h"
#include "ptm_neighbour_ordering.h"
#include "ptm_normalize_vertices.h"


namespace ptm {

typedef struct
{
	double area;
	double dist;
	int index;
	int inner;
	int32_t number;
	double offset[3];
} sorthelper_t;

typedef struct
{
	size_t index;
	int32_t number;
	double area;
	double offset[3];

} solidnbr_t;

static bool sorthelper_compare(sorthelper_t const& a, sorthelper_t const& b)
{
	if (a.area > b.area)
		return true;

	if (a.area < b.area)
		return false;

	if (a.dist < b.dist)
		return true;

	return false;
}

//todo: change voronoi code to return errors rather than exiting
static int calculate_voronoi_face_areas(int num_points, const double (*_points)[3], double* normsq, double max_norm, ptm_voro::voronoicell_neighbor* v, std::vector<int>& nbr_indices, std::vector<double>& face_areas)
{
	const double k = 10 * max_norm;
	v->init(-k,k,-k,k,-k,k);

	for (int i=1;i<num_points;i++)
	{
		double x = _points[i][0] - _points[0][0];
		double y = _points[i][1] - _points[0][1];
		double z = _points[i][2] - _points[0][2];
		v->nplane(x,y,z,normsq[i],i);
	}

	v->neighbors(nbr_indices);
	v->face_areas(face_areas);
	return 0;
}

static int _calculate_neighbour_ordering(void* _voronoi_handle, int num_points, double (*_points)[3], sorthelper_t* data)
{
	assert(num_points <= PTM_MAX_INPUT_POINTS);

	ptm_voro::voronoicell_neighbor* voronoi_handle = (ptm_voro::voronoicell_neighbor*)_voronoi_handle;

	double max_norm = 0;
	double points[PTM_MAX_INPUT_POINTS][3];
	double normsq[PTM_MAX_INPUT_POINTS];

	for (int i=0;i<num_points;i++)
	{
		double x = _points[i][0] - _points[0][0];
		double y = _points[i][1] - _points[0][1];
		double z = _points[i][2] - _points[0][2];
		points[i][0] = x;
		points[i][1] = y;
		points[i][2] = z;

		normsq[i] = x*x + y*y + z*z;
		max_norm = std::max(max_norm, normsq[i]);
	}

	max_norm = sqrt(max_norm);

	std::vector<int> nbr_indices(num_points + 6);
	std::vector<double> face_areas(num_points + 6);
	int ret = calculate_voronoi_face_areas(num_points, points, normsq, max_norm, voronoi_handle, nbr_indices, face_areas);
	if (ret != 0)
		return ret;

	double areas[PTM_MAX_INPUT_POINTS];
	memset(areas, 0, num_points * sizeof(double));
	areas[0] = INFINITY;
	for (size_t i=0;i<nbr_indices.size();i++)
	{
		int index = nbr_indices[i];
		if (index > 0)
			areas[index] = face_areas[i];
	}

	for (int i=0;i<num_points;i++)
	{
		assert(areas[i] == areas[i]);
		data[i].area = areas[i];
		data[i].dist = normsq[i];
		data[i].index = i;
		memcpy(data[i].offset, points[i], 3 * sizeof(double));
	}

	std::sort(data, data + num_points, &sorthelper_compare);
	return ret;
}

int calculate_neighbour_ordering(	void* _voronoi_handle, size_t atom_index, int min_points, int (get_neighbours)(void* vdata, size_t atom_index, int num, size_t* nbr_indices, int32_t* numbers, double (*nbr_pos)[3]), void* nbrlist,
					size_t* indices, double (*points)[3], int32_t* numbers)
{
	size_t nbr_indices[PTM_MAX_INPUT_POINTS];
	int32_t nbr_numbers[PTM_MAX_INPUT_POINTS];
	double nbr_pos[PTM_MAX_INPUT_POINTS][3];
	int num_points = get_neighbours(nbrlist, atom_index, PTM_MAX_INPUT_POINTS, nbr_indices, nbr_numbers, nbr_pos);
	if (num_points < min_points)
		return -1;

	sorthelper_t data[PTM_MAX_INPUT_POINTS];
	int ret = _calculate_neighbour_ordering(_voronoi_handle, num_points, nbr_pos, data);
	if (ret != 0)
		return ret;

	for (int i=0;i<num_points;i++)
	{
		int index = data[i].index;
		indices[i] = nbr_indices[index];
		numbers[i] = nbr_numbers[index];
		memcpy(&points[i], nbr_pos[index], 3 * sizeof(double));
	}

	return num_points;
}

static int find_diamond_neighbours(void* _voronoi_handle, int num_points, double (*_points)[3], size_t* nbr_indices, int32_t* nbr_numbers, int num_solid_nbrs, solidnbr_t* nbrlist)
{
	sorthelper_t data[PTM_MAX_INPUT_POINTS];
	int ret = _calculate_neighbour_ordering(_voronoi_handle, num_points, _points, data);
	if (ret != 0)
		return ret;

	int n = std::min(num_solid_nbrs, num_points - 1);
	for (int i=0;i<n;i++)
	{
		nbrlist[i].index = nbr_indices[data[i+1].index];
		nbrlist[i].area = data[i+1].area;
		nbrlist[i].number = nbr_numbers[data[i+1].index];
		memcpy(nbrlist[i].offset, data[i+1].offset, 3 * sizeof(double));
	}

	return ret;
}

void* voronoi_initialize_local()
{
	ptm_voro::voronoicell_neighbor* ptr = new ptm_voro::voronoicell_neighbor;
	return (void*)ptr;
}

void voronoi_uninitialize_local(void* _ptr)
{
	ptm_voro::voronoicell_neighbor* ptr = (ptm_voro::voronoicell_neighbor*)_ptr;
	delete ptr;
}

int calculate_diamond_neighbour_ordering(	void* _voronoi_handle, size_t atom_index, int (get_neighbours)(void* vdata, size_t atom_index, int num, size_t* nbr_indices, int32_t* numbers, double (*nbr_pos)[3]), void* nbrlist,
						size_t* nbr_indices, double (*points)[3], int32_t* numbers)
{
	const int max_snbrs = 12;
	size_t central_nbr_indices[max_snbrs + 1];
	int32_t central_nbr_numbers[max_snbrs + 1];
	double central_nbr_pos[max_snbrs + 1][3];
	int num_points = get_neighbours(nbrlist, atom_index, max_snbrs + 1, central_nbr_indices, central_nbr_numbers, central_nbr_pos);
	if (num_points < 4)
		return -1;

	solidnbr_t central_solid[max_snbrs];
	int ret = find_diamond_neighbours(_voronoi_handle, num_points, central_nbr_pos, central_nbr_indices, central_nbr_numbers, max_snbrs, central_solid);
	if (ret != 0)
		return ret;


	sorthelper_t data[4 * 6];
	nbr_indices[0] = atom_index;
	numbers[0] = central_nbr_numbers[0];
	memset(&points[0], 0, 3 * sizeof(double));

	std::unordered_set<size_t> claimed;
	claimed.insert(atom_index);
	for (int i=0;i<4;i++)
	{
		size_t inner_index = central_solid[i].index;
		claimed.insert(inner_index);
	}

	int index = 0;
	for (int i=0;i<4;i++)
	{
		size_t inner_index = central_solid[i].index;
		nbr_indices[i+1] = inner_index;
		numbers[i+1] = central_nbr_numbers[i+1];
		memcpy(&points[i+1], central_solid[i].offset, 3 * sizeof(double));

		//----------------------------------------------
		size_t inner_nbr_indices[max_snbrs + 1];
		int32_t inner_nbr_numbers[max_snbrs + 1];
		double inner_nbr_pos[max_snbrs + 1][3];
		int num_points = get_neighbours(nbrlist, inner_index, max_snbrs + 1, inner_nbr_indices, inner_nbr_numbers, inner_nbr_pos);
		if (num_points < 4)
			return -1;

		solidnbr_t inner_solid[max_snbrs];
		ret = find_diamond_neighbours(_voronoi_handle, num_points, inner_nbr_pos, inner_nbr_indices, inner_nbr_numbers, max_snbrs, inner_solid);
		if (ret != 0)
			return ret;
		//----------------------------------------------

		int n = std::min(6, num_points);
		for (int j=0;j<n;j++)
		{
			bool already_claimed = claimed.find(inner_solid[j].index) != claimed.end();
			if (already_claimed)
				continue;

			data[index].inner = i;
			data[index].index = inner_solid[j].index;
			data[index].area = inner_solid[j].area;
			data[index].dist = 0;
			data[index].number = inner_solid[j].number;

			memcpy(data[index].offset, inner_solid[j].offset, 3 * sizeof(double));
			for (int k=0;k<3;k++)
				data[index].offset[k] += central_solid[i].offset[k];

			index++;
		}
	}

	int n = index;
	std::sort(data, data + n, &sorthelper_compare);
	
	int num_found = 0;
	int counts[4] = {0};
	for (int i=0;i<n;i++)
	{
		int inner = data[i].inner;
		int nbr_atom_index = data[i].index;

		bool already_claimed = claimed.find(nbr_atom_index) != claimed.end();
		if (counts[inner] >= 3 || already_claimed)
			continue;

		nbr_indices[1 + 4 + 3 * inner + counts[inner]] = nbr_atom_index;
		numbers[1 + 4 + 3 * inner + counts[inner]] = data[i].number;
		memcpy(points[1 + 4 + 3 * inner + counts[inner]], &data[i].offset, 3 * sizeof(double));
		claimed.insert(nbr_atom_index);

		counts[inner]++;
		num_found++;
		if (num_found >= 12)
			break;
	}

	if (num_found != 12)
		return -1;

	return 0;
}

}

