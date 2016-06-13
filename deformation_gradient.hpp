#ifndef DEFORMATION_GRADIENT_HPP
#define DEFORMATION_GRADIENT_HPP

#include <cstdint>

void calculate_deformation_gradient(int num_points, const double (*ideal_points)[3], int8_t* mapping, double (*normalized)[3], const double (*penrose)[3], double* F, double* res);

extern const double penrose_sc[7][3];
extern const double penrose_fcc[13][3];
extern const double penrose_hcp[13][3];
extern const double penrose_ico[13][3];
extern const double penrose_bcc[15][3];

#endif

