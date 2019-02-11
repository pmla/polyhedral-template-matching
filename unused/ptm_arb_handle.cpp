#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include "ptm_arb_handle.h"


namespace ptm {

void handle_uninit(handle_t** p_handle)
{
	if (p_handle == NULL)
		return;

	handle_t* handle = *p_handle;
	if (handle == NULL)
		return;

	pqueue_free(&handle->pq);
	free(handle->data);
	free(handle->pscratch);
	free(handle->qscratch);
	free(handle);
	*p_handle = NULL;
}

handle_t* handle_init(int num_points, size_t max_memory)
{
	handle_t* handle = (handle_t*)calloc(1, sizeof(handle_t));
	if (handle == NULL)
		return NULL;

	handle->num_points = num_points;
	handle->size = 0;
	handle->capacity = 512;
	handle->max_memory = max_memory;

	handle->pq = pqueue_init();
	handle->data = (node_t*)malloc(sizeof(node_t) * handle->capacity);

	handle->pscratch = (uint8_t*)malloc(sizeof(uint8_t) * num_points * handle->capacity);
	handle->qscratch = (uint8_t*)malloc(sizeof(uint8_t) * num_points * handle->capacity);

	if (handle->pq == NULL || handle->data == NULL || handle->pscratch == NULL || handle->qscratch == NULL)
		handle_uninit(&handle);

	return handle;
}

int handle_add(handle_t* handle, int level, double rmsd, uint8_t* perm_P, uint8_t* perm_Q, double G1, double G2, double* A, bool branch_on_Q, int branch_index)
{
	assert(branch_on_Q == true);

	if (handle->size == handle->capacity)
	{
		if (sizeof(uint8_t) * (handle->capacity * 2) * handle->num_points > handle->max_memory)
			return -2;

		node_t* tmp = NULL;
		if (!(tmp = (node_t*)realloc(handle->data, sizeof(node_t) * handle->capacity * 2)))
			return -1;

		handle->data = tmp;

		uint8_t* temp = NULL;
		if (!(temp = (uint8_t*)realloc(handle->pscratch, handle->num_points*sizeof(uint8_t) * handle->capacity * 2)))
			return -1;

		handle->pscratch = temp;

		if (!(temp = (uint8_t*)realloc(handle->qscratch, handle->num_points*sizeof(uint8_t) * handle->capacity * 2)))
			return -1;

		handle->qscratch = temp;

		handle->capacity *= 2;
	}

	node_t* nd = &handle->data[handle->size];

	nd->level = level;
	nd->rmsd = rmsd;
	nd->G1 = G1;
	nd->G2 = G2;
	memcpy(nd->A, A, 9 * sizeof(double));
	nd->branch_on_Q = branch_on_Q;

	memcpy(&handle->pscratch[handle->size * handle->num_points], perm_P, handle->num_points * sizeof(uint8_t));
	memcpy(&handle->qscratch[handle->size * handle->num_points], perm_Q, handle->num_points * sizeof(uint8_t));
	nd->permutation_index = handle->size * handle->num_points;

	if (branch_on_Q)
	{
		uint8_t* t = &handle->qscratch[handle->size * handle->num_points];
		int temp = t[level];
		t[level] = t[branch_index];
		t[branch_index] = temp;
	}
	else
		assert(1==2);

	int ret = pqueue_insert(handle->pq, rmsd, handle->size);
	if (ret != 0)
		return -1;

	handle->size++;
	return 0;
}

}

