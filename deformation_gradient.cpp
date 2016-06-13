#include <cstdint>
#include "deformation_gradient.hpp"


//sc
#define k_sc 0.5
const double penrose_sc[7][3] = {	
					{0, 0, 0},
					{0, 0, -k_sc},
					{0, 0, k_sc},
					{0, -k_sc, 0},
					{0, k_sc, 0},
					{-k_sc, 0, 0},
					{k_sc, 0, 0},
				};

//fcc
#define k_fcc 0.17677669529663678216
const double penrose_fcc[13][3] = {
					{0, 0, 0},
					{0, k_fcc, k_fcc},
					{0, -k_fcc, -k_fcc},
					{0, k_fcc, -k_fcc},
					{0, -k_fcc, k_fcc},
					{k_fcc, 0, k_fcc},
					{-k_fcc, 0, -k_fcc},
					{k_fcc, 0, -k_fcc},
					{-k_fcc, 0, k_fcc},
					{k_fcc, k_fcc, -0},
					{-k_fcc, -k_fcc, 0},
					{k_fcc, -k_fcc, 0},
					{-k_fcc, k_fcc, -0},
				};

//hcp
#define k_hcp 0.17677669529663678216
const double penrose_hcp[13][3] = {
					{0, 0, 0},
					{k_hcp, 0, k_hcp},
					{-k_hcp/3, -4*k_hcp/3, -k_hcp/3},
					{k_hcp, k_hcp, 0},
					{-k_hcp/3, -k_hcp/3, -4*k_hcp/3},
					{0, k_hcp, k_hcp},
					{-4*k_hcp/3, -k_hcp/3, -k_hcp/3},
					{-k_hcp, k_hcp, -0},
					{0, k_hcp, -k_hcp},
					{k_hcp, 0, -k_hcp},
					{k_hcp, -k_hcp, 0},
					{-k_hcp, 0, k_hcp},
					{0, -k_hcp, k_hcp},
				};

//ico
#define k_ico 0.13143277802974323576
#define phi 1.61803398874989490253
//((1.0 + sqrt(5)) / 2)
const double penrose_ico[13][3] = {
					{0, 0, 0},
					{0, k_ico, phi*k_ico},
					{0, -k_ico, -phi*k_ico},
					{0, k_ico, -phi*k_ico},
					{0, -k_ico, phi*k_ico},
					{-k_ico, -phi*k_ico, -0},
					{k_ico, phi*k_ico, 0},
					{k_ico, -phi*k_ico, 0},
					{-k_ico, phi*k_ico, -0},
					{-phi*k_ico, 0, -k_ico},
					{phi*k_ico, 0, k_ico},
					{phi*k_ico, 0, -k_ico},
					{-phi*k_ico, 0, k_ico},
				};

//bcc
#define k_bcc 0.11543038598460284017
const double penrose_bcc[15][3] = {
					{0, 0, 0},
					{-k_bcc, -k_bcc, -k_bcc},
					{k_bcc, k_bcc, k_bcc},
					{k_bcc, -k_bcc, -k_bcc},
					{-k_bcc, k_bcc, k_bcc},
					{-k_bcc, k_bcc, -k_bcc},
					{k_bcc, -k_bcc, k_bcc},
					{-k_bcc, -k_bcc, k_bcc},
					{k_bcc, k_bcc, -k_bcc},
					{0, 0, -2*k_bcc},
					{0, 0, 2*k_bcc},
					{0, -2*k_bcc, 0},
					{0, 2*k_bcc, 0},
					{-2*k_bcc, 0, 0},
					{2*k_bcc, 0, -0},
				};

void calculate_deformation_gradient(int num_points, const double (*ideal_points)[3], int8_t* mapping, double (*normalized)[3], const double (*penrose)[3], double* F, double* res)
{
	for (int i = 0;i<3;i++)
	{
		for (int j = 0;j<3;j++)
		{
			double acc = 0.0;
			for (int k = 0;k<num_points;k++)
				acc += penrose[k][j] * normalized[mapping[k]][i];

			F[i*3 + j] = acc;
		}
	}

	res[0] = 0;
	res[1] = 0;
	res[2] = 0;

	for (int k = 0;k<num_points;k++)
	{
		for (int i = 0;i<3;i++)
		{
			double acc = 0.0;
			for (int j = 0;j<3;j++)
			{
				acc += F[i*3 + j] * ideal_points[k][j];
			}

			double delta = acc - normalized[mapping[k]][i];
			res[i] += delta * delta;
		}
	}
}

