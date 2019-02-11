#ifndef PTM_ARB_PQUEUE_H
#define PTM_ARB_PQUEUE_H
#include <stdio.h>
#include <string.h>


namespace ptm {

typedef struct
{
	double priority;
	int64_t data;
} pqnode_t;

typedef struct
{
	int size;
	int capacity;
	pqnode_t* data;
} pqueue_t;


pqueue_t* pqueue_init();
void pqueue_free(pqueue_t** p_q);
size_t pqueue_size(pqueue_t* q);
int pqueue_insert(pqueue_t* q, double priority, int64_t data);
int pqueue_pop(pqueue_t* q, double* p_priority, int64_t* p_data);

}

#endif

