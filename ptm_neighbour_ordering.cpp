/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

//todo: normalize vertices

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <set>
#include "ptm_constants.h"
#include "ptm_voronoi_cell.h"
#include "ptm_neighbour_ordering.h"
#include "ptm_normalize_vertices.h"
#include "ptm_linear_assignment.h"


namespace ptm {

typedef struct
{
	double area;
	double dist;
	int index;
	int32_t number;
} sorthelper_t;

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

static double dot_product(double* a, double* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross_product(double* a, double* b, double* c)
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

static double calculate_solid_angle(double* R1, double* R2, double* R3)	//norms of R1-R3 must be 1
{
	double R2R3[3];
	cross_product(R2, R3, R2R3);
	double numerator = dot_product(R1, R2R3);

	double r1r2 = dot_product(R1, R2);
	double r2r3 = dot_product(R2, R3);
	double r3r1 = dot_product(R3, R1);

	double denominator = 1 + r1r2 + r3r1 + r2r3;
	return fabs(2 * atan2(numerator, denominator));
}

//todo: change voronoi code to return errors rather than exiting
static int calculate_voronoi_face_areas(int num_points, const double (*_points)[3], double* normsq, double max_norm, ptm_voro::voronoicell_neighbor* v, bool calc_solid_angles,
						std::vector<int>& nbr_indices, std::vector<double>& face_areas)
{
	const double k = 10 * max_norm;
	v->init(-k,k,-k,k,-k,k);

	for (int i=0;i<num_points;i++)
	{
		double x = _points[i][0];
		double y = _points[i][1];
		double z = _points[i][2];
		v->nplane(x,y,z,normsq[i],i);
	}

	v->neighbors(nbr_indices);

//v->face_areas(face_areas);
	if (!calc_solid_angles)
	{
		v->face_areas(face_areas);
		return 0;
	}
	else
	{
		std::vector<int> face_vertices;
		std::vector<double> vertices;

		v->face_vertices(face_vertices);
		v->vertices(0, 0, 0, vertices);

		size_t num_vertices = vertices.size() / 3;
		for (size_t i=0;i<num_vertices;i++)
		{
			double x = vertices[i * 3 + 0];
			double y = vertices[i * 3 + 1];
			double z = vertices[i * 3 + 2];

			double s = sqrt(x*x + y*y + z*z);
			vertices[i * 3 + 0] /= s;
			vertices[i * 3 + 1] /= s;
			vertices[i * 3 + 2] /= s;
		}

		int num_faces = v->number_of_faces();

#ifdef DEBUG
		printf("number of voronoi faces: %d\n", num_faces);
#endif

//std::vector<double> solids(face_areas.size()+1);

		size_t c = 0;
		for (int current_face=0;current_face<num_faces;current_face++)
		{
			int num = face_vertices[c++];

			int point_index = nbr_indices[current_face];
			if (point_index >= 0)
			{
				double solid_angle = 0;
				int u = face_vertices[c];
				int v = face_vertices[c+1];
				for (int i=2;i<num;i++)
				{
					int w = face_vertices[c+i];
					double omega = calculate_solid_angle(&vertices[u*3], &vertices[v*3], &vertices[w*3]);
					solid_angle += omega;

					v = w;
				}

				face_areas[current_face] = solid_angle;
//solids[current_face] = solid_angle;
				//face_areas[point_index] = solid_angle;
			}

			c += num;
		}

#ifdef DEBUG
printf("\n");
for (int i=0;i<solids.size();i++)
{
	printf("%d\t%f\t%f\n", i, solids[i], face_areas[i]);
}
#endif

		assert(c == face_vertices.size());
		return 0;
	}
}

static int _calculate_neighbour_ordering(void* _voronoi_handle, int num_selected, int* selected, double (*_points)[3], int central, sorthelper_t* data)
{
	assert(num_selected <= PTM_MAX_INPUT_POINTS);

	ptm_voro::voronoicell_neighbor* voronoi_handle = (ptm_voro::voronoicell_neighbor*)_voronoi_handle;

	double max_norm = 0;
	double points[PTM_MAX_INPUT_POINTS][3];
	double normsq[PTM_MAX_INPUT_POINTS];

	for (int i=0;i<num_selected;i++)
	{
		double x = _points[selected[i]][0] - _points[central][0];
		double y = _points[selected[i]][1] - _points[central][1];
		double z = _points[selected[i]][2] - _points[central][2];
		points[i][0] = x;
		points[i][1] = y;
		points[i][2] = z;

		normsq[i] = x*x + y*y + z*z;
		max_norm = std::max(max_norm, normsq[i]);
	}

	max_norm = sqrt(max_norm);

	std::vector<int> nbr_indices(num_selected + 6);
	std::vector<double> face_areas(num_selected + 6);
	int ret = calculate_voronoi_face_areas(num_selected, points, normsq, max_norm, voronoi_handle, true, nbr_indices, face_areas);
	if (ret != 0)
		return ret;

	double areas[PTM_MAX_INPUT_POINTS] = {0};
	for (size_t i=0;i<nbr_indices.size();i++)
	{
		int index = nbr_indices[i];
		if (index >= 0)
			areas[index] = face_areas[i];
	}

	for (int i=0;i<num_selected;i++)
	{
		assert(areas[i] == areas[i]);
		data[i].area = areas[i];
		data[i].dist = normsq[i];
		data[i].index = i;
	}

	std::sort(data, data + num_selected, &sorthelper_compare);
	return ret;
}

int calculate_neighbour_ordering(	void* _voronoi_handle, int num_input_points, double (*input_points)[3], int32_t* input_numbers, ptm::atomicenv_t* output)
{
	sorthelper_t data[PTM_MAX_INPUT_POINTS];

	int num_selected = std::min(18, num_input_points - 1);
	int selected[PTM_MAX_INPUT_POINTS];
	for (int i=0;i<num_selected;i++)
		selected[i] = i+1;

	int ret = _calculate_neighbour_ordering(_voronoi_handle, num_selected, selected, input_points, 0, data);
	if (ret != 0)
		return ret;

	memcpy(output->points[0], input_points[0], 3 * sizeof(double));
	output->ordering[0] = 0;
	output->numbers[0] = input_numbers[0];

	for (int i=0;i<num_selected;i++)
	{
		int index = selected[data[i].index];
		output->ordering[1+i] = index;
		output->numbers[1+i] = input_numbers[index];
		memcpy(output->points[1+i], input_points[index], 3 * sizeof(double));
	}

	for (int i=1+num_selected;i<num_input_points;i++)
	{
		output->ordering[i] = i;
		output->numbers[i] = input_numbers[i];
		memcpy(output->points[i], input_points[i], 3 * sizeof(double));
	}

	return num_input_points;
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

static int euclidean_argsort(int num_input_points, double (*points)[3], int central, int num_inner, int* output_ordering)
{
	sorthelper_t data[PTM_MAX_INPUT_POINTS];

	int num_found = 0;
	for (int i=0;i<num_input_points;i++)
	{
		if (i < num_inner + 1)	//skip
		{
			data[i].area = 0;
			data[i].dist = INFINITY;
		}
		else
		{
			double dx = points[i][0] - points[central][0];
			double dy = points[i][1] - points[central][1];
			double dz = points[i][2] - points[central][2];

			double d = dx*dx + dy*dy + dz*dz;
			data[i].area = 0;
			data[i].dist = d;
			data[i].index = i;
			num_found++;
		}
	}

	std::sort(data, data + num_input_points, &sorthelper_compare);
	for (int i=0;i<num_input_points;i++)
		output_ordering[i] = data[i].index;

	return num_found;
}

#define MAX_INNER 4

int calculate_two_shell_neighbour_ordering(	void* _voronoi_handle, int num_input_points, int num_inner, int num_outer,
						ptm::atomicenv_t* input, ptm::atomicenv_t* output)
{
	assert(num_inner <= MAX_INNER);

//TODO: use rectangular assignment
	//int n = num_inner * num_outer;
int nc = num_input_points - num_inner - 1;

	std::vector< std::vector< double > > cost(nc);
	for (int i=0;i<nc;i++)
		cost[i] = std::vector< double >(nc, 0);

	for (int i=0;i<num_inner;i++)
	{
		int inner_index = 1+i;

		int selected[PTM_MAX_INPUT_POINTS];
		int num_found = euclidean_argsort(num_input_points, input->points, inner_index, num_inner, selected);
		int num_selected = std::min(num_found, 6);

		sorthelper_t data[PTM_MAX_INPUT_POINTS];
		int ret = _calculate_neighbour_ordering(_voronoi_handle, num_selected, selected, input->points, inner_index, data);
		if (ret != 0)
			return ret;

		for (int j=0;j<num_selected;j++)
		{
			int index = selected[data[j].index] - num_inner - 1;
			double area = data[j].area;

			for (int k=0;k<num_outer;k++)
				cost[i * num_outer + k][index] = -area;
		}
	}

	std::vector< int > lmate(nc);
	std::vector< int > rmate(nc);
	ptm::linear_assignment_sum(cost, lmate, rmate);

	for (int i=0;i<num_inner+1;i++)
	{
		output->ordering[i] = input->ordering[i];
		output->numbers[i] = input->numbers[i];
		memcpy(output->points[i], input->points[i], 3 * sizeof(double));
	}

	for (int i=0;i<num_inner;i++)
	{
		for (int j=0;j<num_outer;j++)
		{
			int index = lmate[i * num_outer + j] + 1 + num_inner;

			output->ordering[1 + num_inner + i * num_outer + j] = input->ordering[index];
			output->numbers[1 + num_inner + i * num_outer + j] = input->numbers[index];
			memcpy(output->points[1 + num_inner + i * num_outer + j], input->points[index], 3 * sizeof(double));
		}
	}

	return 0;
}

}

