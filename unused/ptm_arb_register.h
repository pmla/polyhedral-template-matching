#ifndef PTM_ARB_REGISTER_H
#define PTM_ARB_REGISTER_H

namespace ptm {

#define PTM_ARB_MAX_POINTS 25

int register_points(int num_points, double (*P)[3], double (*Q)[3], double max_rmsd, uint8_t* best_permutation, double* p_rmsd, int* p_num_nodes_explored);

}

#endif

