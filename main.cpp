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


//non-diamond data
//#define _MAX_NBRS 24

//diamond data
#define _MAX_NBRS 50

//#define DEBUG

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

static bool get_diamond_neighbours(double (*positions)[3], int32_t* nbrs, int max_nbrs, int central, double (*nbr)[3])
{
	assert(max_nbrs >= 4);

	memcpy(nbr[0], positions[central], 3 * sizeof(double));
	std::set<int> hit;
	hit.insert(central);

int yep[17] = {0};
memset(yep, -1, 17 * sizeof(int));
yep[0] = central;

	std::map< int, int> counts;
	std::map< int, int> pos;
	vector<tuple< double, int, int >> data;
	for (int i=0;i<4;i++)
	{
		int inner = nbrs[central * _MAX_NBRS + i];
yep[1+i] = inner;

		memcpy(nbr[i + 1], positions[inner], 3 * sizeof(double));
		hit.insert(inner);
		counts[inner] = 0;
		pos[inner] = i;

		double* p = positions[inner];
		for (int j=0;j<max_nbrs;j++)
		{
			int outer = nbrs[inner * _MAX_NBRS + j];
			double* q = positions[outer];

			double dx = p[0] - q[0];
			double dy = p[1] - q[1];
			double dz = p[2] - q[2];
			double dist = dx*dx + dy*dy + dz*dz;
			data.push_back(make_tuple(dist, inner, outer));
		}
	}

	sort(data.begin(), data.end(), [](	const tuple<double,int,int>& a,
						const tuple<double,int,int>& b) -> bool
						{
							return std::get<0>(a) < std::get<0>(b);
						});
#ifdef DEBUG
printf("@index: %d\n", central);
#endif
	int num_found = 0;
	for ( auto t = data.begin(); t != data.end(); t++ )
	{
		int inner = std::get<1>(*t);
		int outer = std::get<2>(*t);
		if (counts[inner] >= 3)
			continue;
#ifdef DEBUG
printf("@: %f\t(%d, %d)\t%d\n", std::get<0>(*t), pos[inner], std::get<1>(*t), std::get<2>(*t));
#endif

		if (hit.find(outer) != hit.end())
			continue;

yep[1 + 4 + 3 * pos[inner] + counts[inner]] = outer;

#ifdef DEBUG
printf("@: ");
for (int j=0;j<17;j++)
	if (yep[j] == -1) printf(". ");
	else	printf("%d ", yep[j]);
printf("\n@\n");
#endif

		memcpy(nbr[1 + 4 + 3 * pos[inner] + counts[inner]], positions[outer], 3 * sizeof(double));

		counts[inner]++;
		hit.insert(outer);
		num_found++;

		if (num_found >= 12)
			break;
	}

for (int i=0;i<17;i++)
	assert( yep[i] != -1 );

	return num_found == 12;
}

int main()
{
	ptm_initialize_global();
//uint64_t res = run_tests();
//assert(res == 0);
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
	const int max_nbrs = 18;
	const int num_diamond_nbrs = 16;
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

		double diamond_nbr[num_diamond_nbrs+1][3];
		bool ok = get_diamond_neighbours((double (*)[3])positions, nbrs, max_nbrs, i, diamond_nbr);
		if (!ok)
		{
			//toggle diamond flags;
		}

		int8_t mapping[_MAX_NBRS];
		int32_t type, alloy_type;
		double scale, rmsd, interatomic_distance, lattice_constant;
		double q[4], F[9], F_res[3], U[9], P[9];
		ptm_index(	local_handle, PTM_CHECK_ALL, max_nbrs + 1, nbr, NULL, topological_ordering,
				num_diamond_nbrs + 1, diamond_nbr, NULL,
				&type, &alloy_type, &scale, &rmsd, q, F, F_res, U, P, mapping, &interatomic_distance, &lattice_constant);

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

/*
#ifdef DEBUG
if (type != PTM_MATCH_DCUB || rmsd > 0.05)
	continue;

printf("@type: %d\n", type);
printf("@rmsd: %f\n", rmsd);
for (int j=0;j<17;j++)
{
	printf("!!!%f %f %f\n", nbr[j][0], nbr[j][1], nbr[j][2]);
}
#endif
*/	}

	printf("rmsd sum: %f\n", rmsd_sum);

	free(types);
	free(rmsds);
	free(quats);
	free(positions);
	free(nbrs);
	ptm_uninitialize_local(local_handle);
	return 0;
}

