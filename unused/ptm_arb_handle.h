#ifndef PTM_ARB_HANDLE_H
#define PTM_ARB_HANDLE_H

#include <stdint.h>
#include "ptm_arb_pqueue.h"


namespace ptm {

typedef struct
{
	int level;
	double rmsd;
	double G1;
	double G2;
	double A[9];
	int64_t permutation_index;
	bool branch_on_Q;
	int branch_index;
} node_t;

typedef struct
{
	int64_t size;
	int64_t capacity;
	size_t max_memory;
	int64_t num_points;
	pqueue_t* pq;
	node_t* data;
	uint8_t* pscratch;
	uint8_t* qscratch;
} handle_t;

void handle_uninit(handle_t** p_handle);
handle_t* handle_init(int num_points, size_t max_memory);
int handle_add(handle_t* handle, int level, double rmsd, uint8_t* perm_P, uint8_t* perm_Q, double G1, double G2, double* A, bool branch_on_Q, int branch_index);

}

#endif

