#include <stdint.h>
#include <stdbool.h>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cassert>
#include <algorithm>
#include "ptm_functions.h"
#include "unittest.hpp"

using namespace std;

//counts: [0 10726 66758 16207 668 339181 569499 42834 285213 ]
//counts: [0 10491 65040 15990 660 324483 602098 40754 271570 ]

//diamond test results:
//counts: [38 7433 45410 8685 479 261154 657373 68510 282004 ]
//rmsd sum: 302846.681637


#define _MAX_NBRS 50	//diamond
//#define _MAX_NBRS 24	//fcc, other
//#define _MAX_NBRS 6	//graphene


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

int main()
{
	ptm_initialize_global();
	uint64_t res = ptm::run_tests();
	assert(res == 0);
	//printf("=========================================================\n");
	//printf("unit test result: %lu\n", res);
	//return 0;

	size_t fsize = 0;
	int32_t* nbrs = NULL;
	double (*positions)[3] = NULL;
	//int ret = read_file((char*)"test_data/FeCu_positions.dat", (uint8_t**)&positions, &fsize);
	//int ret = read_file((char*)"test_data/fcc_positions.dat", (uint8_t**)&positions, &fsize);
	int ret = read_file((char*)"test_data/diamond_pos.dat", (uint8_t**)&positions, &fsize);
	//int ret = read_file((char*)"test_data/graphene_pos.dat", (uint8_t**)&positions, &fsize);
	if (ret != 0)
		return -1;

	//ret = read_file((char*)"test_data/FeCu_nbrs.dat", (uint8_t**)&nbrs, &fsize);
	//ret = read_file((char*)"test_data/fcc_nbrs.dat", (uint8_t**)&nbrs, &fsize);
	ret = read_file((char*)"test_data/diamond_nbrs.dat", (uint8_t**)&nbrs, &fsize);
	//ret = read_file((char*)"test_data/graphene_nbrs.dat", (uint8_t**)&nbrs, &fsize);
	if (ret != 0)
		return -1;

	int num_atoms = fsize / (_MAX_NBRS * sizeof(int32_t));

	//assert(num_atoms == 88737);
	printf("num atoms: %d\n", num_atoms);

	int8_t* types = (int8_t*)calloc(sizeof(int8_t), num_atoms);
	double* rmsds = (double*)calloc(sizeof(double), num_atoms);
	double* quats = (double*)calloc(4 * sizeof(double), num_atoms);
	int counts[9] = {0};

	ptm_local_handle_t local_handle = ptm_initialize_local();

	double rmsd_sum = 0.0;
	for (int i=0;i<num_atoms;i++)
	{
		double nbr_pos[PTM_MAX_INPUT_POINTS][3];
		int32_t numbers[PTM_MAX_INPUT_POINTS];

		memcpy(nbr_pos[0], positions[i], 3 * sizeof(double));
		if (numbers != NULL)
			numbers[0] = 0;

		int n = std::min(PTM_MAX_INPUT_POINTS - 1, _MAX_NBRS);
		for (int j=0;j<n;j++)
		{
			size_t index = nbrs[i * _MAX_NBRS + j];
			if (numbers != NULL)
				numbers[j+1] = 0;

			memcpy(nbr_pos[j+1], positions[index], 3 * sizeof(double));
		}

		int8_t output_indices[PTM_MAX_INPUT_POINTS];
		int32_t type, alloy_type;
		double scale, rmsd, interatomic_distance, lattice_constant;
		double q[4], F[9], F_res[3], U[9], P[9];
		const double (*best_template)[3];

		ptm_index(	local_handle, n, nbr_pos, NULL, PTM_CHECK_ALL, true,
				&type, &alloy_type, &scale, &rmsd, q, F, F_res, U, P, &interatomic_distance, &lattice_constant, &best_template, output_indices);
//{
//	printf("#scale %f\n", scale);
//	printf("#quat %f %f %f %f\n", q[0], q[1], q[2], q[3]);
//	for (int j=0;j<32;j++)
//		printf("!!!%f %f %f\n", nbr_pos[j][0], nbr_pos[j][1], nbr_pos[j][2]);
//}

		types[i] = type;
		rmsds[i] = rmsd;
		memcpy(&quats[i*4], q, 4 * sizeof(double));

		counts[type]++;
		if (type != PTM_MATCH_NONE)
			rmsd_sum += rmsd;

		if (i % 1000 == 0 || i == num_atoms - 1 || 0)
		{
			printf("counts: [");
			for (int j = 0;j<9;j++)
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

