#include <cstdlib>
#include <cmath>
#include <cassert>
#include "voronoi/cell.hpp"
using namespace voro;
using namespace std;

#define MAX_POINTS 19


#ifdef __cplusplus
extern "C" {
#endif

int calculate_neighbour_ordering(voronoicell_neighbor* voronoi_handle, int num_points, const double (*_points)[3], int8_t* ordering);

void* ptm_initialize_local();
void ptm_uninitialize_local(voronoicell_neighbor* ptr);


#ifdef __cplusplus
}
#endif


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

//todo: change voronoi code to return errors rather than exiting
static int calculate_voronoi_face_areas(int num_points, const double (*_points)[3], double* normsq, double max_norm, voronoicell_neighbor* v, vector<int>& nbr_indices, vector<double>& face_areas)
{
	const double k = 1000 * max_norm;
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

int calculate_neighbour_ordering(voronoicell_neighbor* voronoi_handle, int num_points, const double (*_points)[3], int8_t* ordering)
{
	assert(num_points <= MAX_POINTS);

	double max_norm = 0;
	double points[num_points][3];
	double normsq[num_points];
	for (int i = 0;i<num_points;i++)
	{
		double x = _points[i][0] - _points[0][0];
		double y = _points[i][1] - _points[0][1];
		double z = _points[i][2] - _points[0][2];
		points[i][0] = x;
		points[i][1] = y;
		points[i][2] = z;

		normsq[i] = x*x + y*y + z*z;
		max_norm = max(max_norm, normsq[i]);
#ifdef DEBUG
		printf("point %d: %f\t%f\t%f\t%f\n", i, x, y, z, x*x + y*y + z*z);
#endif
	}

	max_norm = sqrt(max_norm);

	vector<int> nbr_indices(num_points+6);
	vector<double> face_areas(num_points+6);
	int ret = calculate_voronoi_face_areas(num_points, points, normsq, max_norm, voronoi_handle, nbr_indices, face_areas);
	if (ret != 0)
		return ret;

	double areas[num_points] = {0};
	areas[0] = INFINITY;
	for (size_t i=0;i<nbr_indices.size();i++)
	{
		int index = nbr_indices[i];
		if (index > 0)
			areas[index] = face_areas[i];
	}

	sorthelper_t data[MAX_POINTS];
	for (int i=0;i<num_points;i++)
	{
		assert(areas[i] == areas[i]);
		data[i].area = areas[i];
		data[i].dist = normsq[i];
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

void* ptm_initialize_local()
{
	voronoicell_neighbor* ptr = new voronoicell_neighbor;
	return (void*)ptr;
}

void ptm_uninitialize_local(voronoicell_neighbor* ptr)
{
	delete ptr;
}

