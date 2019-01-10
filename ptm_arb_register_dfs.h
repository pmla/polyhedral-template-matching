#ifndef PTM_ARB_REGISTER_DFS_H
#define PTM_ARB_REGISTER_DFS_H

#include "ptm_arb_treefilter.h"

namespace ptm {

int register_points_dfs(int num_points, double (*P)[3], double (*Q)[3], double max_rmsd, const treefilter_t* filter,
				uint8_t* best_permutation, double* p_rmsd, int* p_num_nodes_explored);

}

#endif

