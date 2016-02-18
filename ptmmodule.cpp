#include <Python.h>
#include <ndarraytypes.h>
#include "index_PTM.h"


#define MIN_NBRS 6

extern "C" {

static PyObject* error(PyObject* type, char* msg)
{
	PyErr_SetString(type, msg);
	return NULL;
}

static PyObject* index_structure(PyObject* self, PyObject* args, PyObject* kw)
{
	PyArrayObject* obj_pos = NULL;
	PyArrayObject* obj_num = NULL;
	PyObject* obj_types = NULL;

	static char* argnames[] = {"rel", "numbers", "structures", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kw, "O|OO", argnames, &obj_pos, &obj_num, &obj_types))
		return NULL;

	int num_nbrs = 15;

	if (PyArray_NDIM(obj_pos) != 2			//two-dimensional
		|| PyArray_DIM(obj_pos, 0) != num_nbrs	//need 14 nearest neighbours + central atom
		|| PyArray_DIM(obj_pos, 1) != 3		//second dim is 3
		|| PyArray_TYPE(obj_pos) != NPY_DOUBLE	//array of double
		|| !PyArray_ISCARRAY_RO(obj_pos))	//contiguous etc.
		return error(PyExc_TypeError, "neighbour array must be 14x3 double array");

	bool check_alloys = obj_num != NULL;
	if (check_alloys)
	{
		//if (PyArray_NDIM(obj_num) != 1			//one-dimensional
		//	|| PyArray_DIM(obj_num, 0) != num_nbrs	//need 14 nearest neighbours + central atom
		//	|| PyArray_TYPE(obj_num) != NPY_INT	//array of ints
		//	|| !PyArray_ISCARRAY_RO(obj_num))	//contiguous etc.
		//	return error(PyExc_TypeError, "numbers array must be 14x1 integer array");

		if (PyArray_NDIM(obj_num) != 1)			//one-dimensional
			return error(PyExc_TypeError, "numbers array must be 1-dimensional");

		if (PyArray_DIM(obj_num, 0) != num_nbrs)	//need 14 nearest neighbours + central atom
			return error(PyExc_TypeError, "numbers array must be contain 15 elements");

		if (PyArray_TYPE(obj_num) != NPY_INT)	//array of ints
			return error(PyExc_TypeError, "numbers array must be have dtype NPY_INT");

		if (!PyArray_ISCARRAY_RO(obj_num))	//contiguous etc.
			return error(PyExc_TypeError, "numbers array must be contiguous array");
	}

	int32_t flags = 0;
	if (obj_types == NULL)
	{
		flags = PTM_CHECK_ALL;
	}
	else
	{
		if (!PyTuple_Check(obj_types))
			return error(PyExc_TypeError, "types must be a tuple of strings");

		int num_types = PyTuple_Size(obj_types);

		int i = 0;
		for (i=0;i<num_types;i++)
		{
			PyObject* obj_type = PyTuple_GetItem(obj_types, i);
			if (obj_type == NULL)
				return NULL;

			if (!PyString_Check(obj_type))
				return error(PyExc_TypeError, "type is not a string");

			char* type = PyBytes_AsString(obj_type);
			if (type == NULL)
				return NULL;

			if (strcmp(type, "sc") == 0)
				flags |= PTM_CHECK_SC;
			else if (strcmp(type, "fcc") == 0)
				flags |= PTM_CHECK_FCC;
			else if (strcmp(type, "hcp") == 0)
				flags |= PTM_CHECK_HCP;
			else if (strcmp(type, "ico") == 0)
				flags |= PTM_CHECK_ICO;
			else if (strcmp(type, "bcc") == 0)
				flags |= PTM_CHECK_BCC;
			else
				return error(PyExc_ValueError, "unrecognized type string");
		}
	}

	double* pos = (double*)PyArray_DATA(obj_pos);
	if (pos == NULL)
		return NULL;

	int32_t* numbers = NULL;
	if (check_alloys)
	{
		numbers = (int32_t*)PyArray_DATA(obj_num);
		if (numbers == NULL)
			return NULL;
	}

	int32_t type, alloy_type;
	double scale, rmsd;
	double q[4], F[9], lstsq_residual[3], U[9], P[9];
	index_PTM(num_nbrs, pos, numbers, flags, &type, &alloy_type, &scale, &rmsd, q, F, lstsq_residual, U, P);
	if (type == PTM_MATCH_NONE)
		return Py_BuildValue("iidd()()()()()", PTM_MATCH_NONE, PTM_ALLOY_NONE, INFINITY, INFINITY);

	return Py_BuildValue("iidd(dddd)(ddddddddd)(ddd)(ddddddddd)(ddddddddd)",
											type, alloy_type, rmsd, scale,
											q[0], q[1], q[2], q[3],
											F[0], F[1], F[2], F[3], F[4], F[5], F[6], F[7], F[8],
											lstsq_residual[0], lstsq_residual[1], lstsq_residual[2],
											P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8],
											U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7], U[8]);


	/*if (res->type == PTM_MATCH_SC)
		return Py_BuildValue("iidd(dddd)(ddddddddd)(ddd)(ddddddddd)(ddddddddd)(iiiiiii)",
											res->type, res->alloy_type, res->rmsd, res->scale,
											q[0], q[1], q[2], q[3],
											F[0], F[1], F[2], F[3], F[4], F[5], F[6], F[7], F[8],
											residual[0], residual[1], residual[2],
											P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8],
											U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7], U[8],
											b[0], b[1], b[2], b[3], b[4], b[5], b[6]);

	else if (res->type == PTM_MATCH_FCC || res->type == PTM_MATCH_HCP || res->type == PTM_MATCH_ICO)
		return Py_BuildValue("iidd(dddd)(ddddddddd)(ddd)(ddddddddd)(ddddddddd)(iiiiiiiiiiiii)",
											res->type, res->alloy_type, res->rmsd, res->scale,
											q[0], q[1], q[2], q[3],
											F[0], F[1], F[2], F[3], F[4], F[5], F[6], F[7], F[8],
											residual[0], residual[1], residual[2],
											P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8],
											U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7], U[8],
											b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12]);
	else
		return Py_BuildValue("iidd(dddd)(ddddddddd)(ddd)(ddddddddd)(ddddddddd)(iiiiiiiiiiiiiii)",
											res->type, res->alloy_type, res->rmsd, res->scale,
											q[0], q[1], q[2], q[3],
											F[0], F[1], F[2], F[3], F[4], F[5], F[6], F[7], F[8],
											residual[0], residual[1], residual[2],
											P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8],
											U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7], U[8],
											b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14]);*/
}

static PyMethodDef PTMModuleMethods[] =
{
	{"index_structure", (PyCFunction)index_structure, METH_VARARGS | METH_KEYWORDS, "determine the structure of the atom"},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initptmmodule(void)
{
	(void) Py_InitModule("ptmmodule", PTMModuleMethods);
	initialize_PTM();
}

}
