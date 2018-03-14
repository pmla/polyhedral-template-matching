#include <cstdint>
#include <cstdbool>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cassert>

#include <map>
#include <set>
#include <tuple>
#include <algorithm>

#include "ptm_functions.h"
#include "unittest.hpp"

using namespace std;

#define ACTIVE_DIAMOND

#ifdef ACTIVE_DIAMOND
	#define _MAX_NBRS 50
#else
	#define _MAX_NBRS 24
#endif

static int read_file(const char* path, uint8_t** p_buf, size_t* p_fsize)
{
	size_t fsize = 0, num_read = 0;
	uint8_t* buf = NULL;

	FILE* fin = fopen(path, "rb");
	if (fin == NULL)
		return -1;

	int ret = fseek(fin, 0, SEEK_END);
	if (ret != 0)
		goto cleanup;

	fsize = ftell(fin);
	ret = fseek(fin, 0, SEEK_SET);
	if (ret != 0)
		goto cleanup;

	buf = (uint8_t*)malloc(fsize);
	if (buf == NULL)
	{
		ret = -1;
		goto cleanup;
	}

	num_read = fread(buf, 1, fsize, fin);
	if (num_read != fsize)
	{
		ret = -1;
		goto cleanup;
	}

cleanup:
	if (ret != 0)
	{
		free(buf);
		buf = NULL;
	}

	*p_fsize = fsize;
	*p_buf = buf;
	fclose(fin);
	return ret;
}

static void get_neighbours(double (*positions)[3], int32_t* nbrs, int max_nbrs, int i, double (*nbr)[3])
{
	memcpy(nbr[0], positions[i], 3 * sizeof(double));

	for (int j=0;j<max_nbrs;j++)
	{
		int index = nbrs[i * _MAX_NBRS + j];
		memcpy(nbr[j+1], positions[index], 3 * sizeof(double));
	}
}

//todo: put all unit tests back in

int main()
{
	ptm_initialize_global();
	uint64_t res = run_tests();
	assert(res == 0);
	//printf("=========================================================\n");
	//printf("unit test result: %lu\n", res);
	//return 0;

	size_t fsize = 0;
	int32_t* nbrs = NULL;
	double* positions = NULL;
	//int ret = read_file((char*)"test_data/FeCu_positions.dat", (uint8_t**)&positions, &fsize);
	//int ret = read_file((char*)"test_data/fcc_positions.dat", (uint8_t**)&positions, &fsize);
	int ret = read_file((char*)"test_data/diamond_pos.dat", (uint8_t**)&positions, &fsize);
	if (ret != 0)
		return -1;

	//ret = read_file((char*)"test_data/FeCu_nbrs.dat", (uint8_t**)&nbrs, &fsize);
	//ret = read_file((char*)"test_data/fcc_nbrs.dat", (uint8_t**)&nbrs, &fsize);
	ret = read_file((char*)"test_data/diamond_nbrs.dat", (uint8_t**)&nbrs, &fsize);
	if (ret != 0)
		return -1;

	int num_atoms = fsize / (_MAX_NBRS * sizeof(int32_t));

#ifdef ACTIVE_DIAMOND
	const int max_nbrs = 34;
#else
	const int max_nbrs = 19;
#endif

	//assert(num_atoms == 88737);
	printf("num atoms: %d\n", num_atoms);

	int8_t* types = (int8_t*)calloc(sizeof(int8_t), num_atoms);
	double* rmsds = (double*)calloc(sizeof(double), num_atoms);
	double* quats = (double*)calloc(4 * sizeof(double), num_atoms);
	int counts[8] = {0};

	ptm_local_handle_t local_handle = ptm_initialize_local();

//num_atoms = std::min(num_atoms, 100000);
	bool topological_ordering = true;
	double rmsd_sum = 0.0;
	for (int i=0;i<num_atoms;i++)
	{
		double nbr[max_nbrs+1][3];
		get_neighbours((double (*)[3])positions, nbrs, max_nbrs, i, nbr);

		int8_t mapping[_MAX_NBRS];
		int32_t type, alloy_type;
		double scale, rmsd, interatomic_distance, lattice_constant;
		double q[4], F[9], F_res[3], U[9], P[9];
		ptm_index(	local_handle, PTM_CHECK_ALL, max_nbrs + 1, nbr, NULL, topological_ordering,
				&type, &alloy_type, &scale, &rmsd, q, F, F_res, U, P, mapping, &interatomic_distance, &lattice_constant);

/*if (type == PTM_MATCH_BCC && rmsd < 0.1)
{
	printf("#scale %f\n", scale);
	printf("#quat %f %f %f %f\n", q[0], q[1], q[2], q[3]);
	for (int j=0;j<15;j++)
		printf("!!!%f %f %f\n", nbr[j][0], nbr[j][1], nbr[j][2]);
}*/

		types[i] = type;
		rmsds[i] = rmsd;
		memcpy(&quats[i*4], q, 4 * sizeof(double));

		counts[type]++;
		if (type != PTM_MATCH_NONE)
			rmsd_sum += rmsd;

		if (i % 1000 == 0 || i == num_atoms - 1 || 0)
		{
			printf("counts: [");
			for (int j = 0;j<8;j++)
				printf("%d ", counts[j]);
			printf("]\n");
		}
	}

	printf("rmsd sum: %f\n", rmsd_sum);

	free(types);
	free(rmsds);
	free(quats);
	free(positions);
	free(nbrs);
	ptm_uninitialize_local(local_handle);
	return 0;
}

