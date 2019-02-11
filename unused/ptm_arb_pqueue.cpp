#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "ptm_arb_pqueue.h"


#define LEFT(i)   (((i+1) * 2) - 1)
#define RIGHT(i)  (((i+1) * 2) + 0)
#define PARENT(i) ((i) >> 1)


namespace ptm {

pqueue_t* pqueue_init()
{
	pqueue_t* q = (pqueue_t*)calloc(1, sizeof(pqueue_t));
	if (q == NULL)
		return NULL;

	q->capacity = 512;
	q->data = (pqnode_t*)malloc(sizeof(pqnode_t) * q->capacity);
	if (q->data == NULL)
	{
		free(q);
		q = NULL;
	}

	return q;
}

void pqueue_free(pqueue_t** p_q)
{
	if (p_q == NULL)
		return;

	pqueue_t* q = *p_q;
	if (q == NULL)
		return;

	free(q->data);
	free(q);
	*p_q = NULL;
}

size_t pqueue_size(pqueue_t* q)
{
	return q->size;
}

int pqueue_insert(pqueue_t* q, double priority, int64_t data)
{
	if (q == NULL)
		return -1;

	if (q->size == q->capacity)
	{
		pqnode_t* tmp;
		if (!(tmp = (pqnode_t*)realloc(q->data, sizeof(pqnode_t) * q->capacity * 2)))
			return -1;

		q->data = tmp;
		q->capacity *= 2;
	}

	pqnode_t nd = {priority, data};

	int index = q->size++;
	while (index > 0)
	{
		int parent = PARENT(index);
		if (nd.priority < q->data[parent].priority)
		{
			memcpy(&q->data[index], &q->data[parent], sizeof(pqnode_t));
			index = parent;
		}
		else
			break;
	}

	memcpy(&q->data[index], &nd, sizeof(pqnode_t));
	return 0;
}

int pqueue_pop(pqueue_t* q, double* p_priority, int64_t* p_data)
{
	if (q == NULL || q->size == 0)
		return -1;

	*p_priority = q->data[0].priority;
	*p_data = q->data[0].data;

	pqnode_t nd;
	memcpy(&nd, &q->data[q->size-1], sizeof(pqnode_t));
	q->size--;

	int index = 0;
	while (index < q->size)
	{
		int left = LEFT(index);
		int right = RIGHT(index);
		if (right < q->size && q->data[right].priority < q->data[left].priority && nd.priority > q->data[right].priority)
		{
			memcpy(&q->data[index], &q->data[right], sizeof(pqnode_t));
			index = right;
		}
		else if (left < q->size && nd.priority > q->data[left].priority)
		{
			memcpy(&q->data[index], &q->data[left], sizeof(pqnode_t));
			index = left;
		}
		else
			break;
	}

	memcpy(&q->data[index], &nd, sizeof(pqnode_t));
	return 0;
}

}

