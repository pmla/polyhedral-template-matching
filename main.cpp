#include <cassert>
#include <cstring>
#include <cfloat>
#include <vector>
#include <algorithm>
#include <fstream>
#include "point.hpp"
#include "index_PTM.h"

#define MAX_NBRS 14

char* read_binary_file(char* filename, std::streamsize& size)
{
	std::ifstream file(filename, std::ios::binary);
	file.seekg(0, std::ios::end);
	size = file.tellg();
	file.seekg(0, std::ios::beg);

	char* buffer = new char[size];
	//if (file.read(buffer.data(), size))
	file.read(buffer, size);
	return buffer;
}

int dump_buffer(const char* path, size_t nbytes, uint8_t* buf)
{
	FILE* fout = fopen(path, "wb");
	if (fout == NULL)
		return -1;

	size_t nwritten = fwrite((char*)buf, 1, nbytes, fout);
	fclose(fout);

	if (nwritten != nbytes)
		return -1;

	return 0;
}

static bool point_comparator(const Point &a, const Point &b)
{
	double sa = a.x*a.x + a.y*a.y + a.z*a.z;
	double sb = b.x*b.x + b.y*b.y + b.z*b.z;
	return sa < sb;
}

static std::vector<Point> get_neighbours(Point* positions, int32_t* nbrs, int max_nbrs, int i)
{
	Point p = positions[i];

	std::vector<Point> nbr;
	nbr.push_back(Point(0, 0, 0));
	for (int j=0;j<max_nbrs;j++)
	{
		int index = nbrs[i * max_nbrs + j];

		Point q = positions[index];
		Point d = Point(q.x - p.x, q.y - p.y, q.z - p.z);
		nbr.push_back(d);
	}

	std::sort(nbr.begin(), nbr.end(), point_comparator);
	return nbr;
}

int main()
{
	std::streamsize size;
	Point* positions = (Point*)read_binary_file((char*)"FeCu_positions.dat", size);
	int32_t* nbrs = (int32_t*)read_binary_file((char*)"FeCu_nbrs.dat", size);

	const int max_nbrs = MAX_NBRS;
	int num_atoms = size / (max_nbrs * sizeof(int32_t));

	initialize_PTM();
	int8_t* types = (int8_t*)calloc(sizeof(int8_t), num_atoms);
	double* rmsds = (double*)calloc(sizeof(double), num_atoms);
	double* quats = (double*)calloc(4 * sizeof(double), num_atoms);
	int counts[6] = {0};

	double rmsd_sum = 0.0;
	for (int i=0;i<num_atoms;i++)
	{
		std::vector<Point> nbr = get_neighbours(positions, nbrs, max_nbrs, i);

		int32_t type, alloy_type;
		double scale, rmsd;
		double q[4], F[9], F_res[3], U[9], P[9];
		index_PTM(nbr.size(), &nbr[0].x, NULL, PTM_CHECK_ALL, &type, &alloy_type, &scale, &rmsd, q, F, F_res, U, P);

		types[i] = type;
		rmsds[i] = rmsd;
		memcpy(&quats[i*4], q, 4 * sizeof(double));

		counts[type]++;
		if (type != PTM_MATCH_NONE)
			rmsd_sum += rmsd;

		if (i % 1000 == 0 || i == num_atoms - 1)
		{
			printf("counts: [");
			for (int j = 0;j<6;j++)
				printf("%d ", counts[j]);
			printf("]\n");
		}
	}

	printf("rmsd sum: %f\n", rmsd_sum);

	/*int ret = dump_buffer("PTM_types.dat", sizeof(int8_t) * num_atoms, (uint8_t*)types);
	ret |= dump_buffer("PTM_rmsds.dat", sizeof(double) * num_atoms, (uint8_t*)rmsds);
	ret |= dump_buffer("PTM_quats.dat", 4 * sizeof(double) * num_atoms, (uint8_t*)quats);
	if (ret != 0)
		printf("couldn't dump results\n");*/

	free(types);
	free(rmsds);
	free(quats);
	delete[] positions;
	delete[] nbrs;
	return 0;
}

