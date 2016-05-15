#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "index_ptm.h"
#include "unittest.h"

#define MAX_NBRS 24


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
		int index = nbrs[i * MAX_NBRS + j];
		memcpy(nbr[j+1], positions[index], 3 * sizeof(double));
	}
}

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
	int ret = read_file((char*)"test_data/FeCu_positions.dat", (uint8_t**)&positions, &fsize);
	//int ret = read_file((char*)"test_data/fcc_positions.dat", (uint8_t**)&positions, &fsize);
	if (ret != 0)
		return -1;

	ret = read_file((char*)"test_data/FeCu_nbrs.dat", (uint8_t**)&nbrs, &fsize);
	//ret = read_file((char*)"test_data/fcc_nbrs.dat", (uint8_t**)&nbrs, &fsize);
	if (ret != 0)
		return -1;

	int num_atoms = fsize / (MAX_NBRS * sizeof(int32_t));
	const int max_nbrs = 18;
	//assert(num_atoms == 88737);
	printf("num atoms: %d\n", num_atoms);

	int8_t* types = (int8_t*)calloc(sizeof(int8_t), num_atoms);
	double* rmsds = (double*)calloc(sizeof(double), num_atoms);
	double* quats = (double*)calloc(4 * sizeof(double), num_atoms);
	int counts[6] = {0};

#ifdef DEBUG
extern int failcount;
#endif

	ptm_local_handle_t local_handle = ptm_initialize_local();

	bool topological_ordering = true;
	double rmsd_sum = 0.0;
	for (int j=0;j<10;j++)
	for (int i=0;i<num_atoms;i++)
	{
		double nbr[max_nbrs+1][3];
		get_neighbours((double (*)[3])positions, nbrs, max_nbrs, i, nbr);

		int8_t mapping[MAX_NBRS];
		int32_t type, alloy_type;
		double scale, rmsd;
		double q[4], F[9], F_res[3], U[9], P[9];
		ptm_index(local_handle, max_nbrs + 1, nbr[0], NULL, PTM_CHECK_ALL, topological_ordering, &type, &alloy_type, &scale, &rmsd, q, F, F_res, U, P, mapping);

		types[i] = type;
		rmsds[i] = rmsd;
		memcpy(&quats[i*4], q, 4 * sizeof(double));

		counts[type]++;
		if (type != PTM_MATCH_NONE)
			rmsd_sum += rmsd;

		if (i % 1000 == 0 || i == num_atoms - 1 || 0)
		{
			printf("counts: [");
			for (int j = 0;j<6;j++)
				printf("%d ", counts[j]);
			printf("]\n");
		}
	}

	printf("rmsd sum: %f\n", rmsd_sum);

#ifdef DEBUG
	printf("failcount: %d\n", failcount);
#endif

	free(types);
	free(rmsds);
	free(quats);
	free(positions);
	free(nbrs);
	ptm_uninitialize_local(local_handle);
	return 0;
}

