#include <Python.h>
#include <ndarraytypes.h>
#include <arrayobject.h>
#include <stdbool.h>
#include "index_PTM.h"
#include "unittest.h"


#define MIN_NBRS 6

#ifdef __cplusplus
extern "C" {
#endif

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
	PyObject* obj_strains = NULL;

	static char* argnames[] = {"rel", "numbers", "structures", "calculate_strains", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kw, "O|OOO", argnames, &obj_pos, &obj_num, &obj_types, &obj_strains))
		return NULL;

	int num_nbrs = 15;

	if (PyArray_NDIM(obj_pos) != 2			//two-dimensional
		|| PyArray_DIM(obj_pos, 0) != num_nbrs	//need 15 nearest neighbours + central atom
		|| PyArray_DIM(obj_pos, 1) != 3		//second dim is 3
		|| PyArray_TYPE(obj_pos) != NPY_DOUBLE	//array of double
		|| !PyArray_ISCARRAY_RO(obj_pos))	//contiguous etc.
		return error(PyExc_TypeError, "neighbour array must be 15x3 double array");

	bool check_alloys = obj_num != NULL;
	if (check_alloys)
	{
		if (PyArray_NDIM(obj_num) != 1)			//one-dimensional
			return error(PyExc_TypeError, "numbers array must be 1-dimensional");

		if (PyArray_DIM(obj_num, 0) != num_nbrs)	//need 14 nearest neighbours + central atom
			return error(PyExc_TypeError, "numbers array must be contain 15 elements");

		if (PyArray_TYPE(obj_num) != NPY_INT)		//array of ints
			return error(PyExc_TypeError, "numbers array must be have dtype NPY_INT");

		if (!PyArray_ISCARRAY_RO(obj_num))		//contiguous etc.
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

	bool calculate_strains = false;
	if (obj_strains != NULL)
	{
		int ret = PyObject_IsTrue(obj_strains);
		if (ret == -1)
			return NULL;

		calculate_strains = ret == 1;
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

	npy_intp dims_3[2] = {3};
	npy_intp dims_4[2] = {4};
	npy_intp dims_3_3[2] = {3, 3};


	if (calculate_strains)
	{
		index_PTM(num_nbrs, pos, numbers, flags, &type, &alloy_type, &scale, &rmsd, q, F, lstsq_residual, U, P);
		if (type == PTM_MATCH_NONE)
			return Py_BuildValue("iiddOOOOO", PTM_MATCH_NONE, PTM_ALLOY_NONE, INFINITY, INFINITY, Py_None, Py_None, Py_None, Py_None, Py_None);

		PyObject* arr_q = PyArray_SimpleNew(1, dims_4, NPY_DOUBLE);
		PyObject* arr_res = PyArray_SimpleNew(1, dims_3, NPY_DOUBLE);
		PyObject* arr_F = PyArray_SimpleNew(2, dims_3_3, NPY_DOUBLE);
		PyObject* arr_P = PyArray_SimpleNew(2, dims_3_3, NPY_DOUBLE);
		PyObject* arr_U = PyArray_SimpleNew(2, dims_3_3, NPY_DOUBLE);

		memcpy(PyArray_DATA((PyArrayObject*)arr_res), lstsq_residual, 3 * sizeof(double));
		memcpy(PyArray_DATA((PyArrayObject*)arr_F), F, 9 * sizeof(double));
		memcpy(PyArray_DATA((PyArrayObject*)arr_P), P, 9 * sizeof(double));
		memcpy(PyArray_DATA((PyArrayObject*)arr_U), U, 9 * sizeof(double));
		memcpy(PyArray_DATA((PyArrayObject*)arr_q), q, 4 * sizeof(double));

		PyObject* result = Py_BuildValue("iiddOOOOO", type, alloy_type, rmsd, scale, arr_q, arr_F, arr_res, arr_P, arr_U);

		Py_DECREF(arr_q);
		Py_DECREF(arr_res);
		Py_DECREF(arr_F);
		Py_DECREF(arr_P);
		Py_DECREF(arr_U);
		return result;
	}
	else
	{
		index_PTM(num_nbrs, pos, numbers, flags, &type, &alloy_type, &scale, &rmsd, q, NULL, NULL, NULL, NULL);
		if (type == PTM_MATCH_NONE)
			return Py_BuildValue("iiddO", PTM_MATCH_NONE, PTM_ALLOY_NONE, INFINITY, INFINITY, Py_None);

		PyObject* arr_q = PyArray_SimpleNew(1, dims_4, NPY_DOUBLE);
		memcpy(PyArray_DATA((PyArrayObject*)arr_q), q, 4 * sizeof(double));
		PyObject* result = Py_BuildValue("iiddO", type, alloy_type, rmsd, scale, arr_q);
		Py_DECREF(arr_q);
		return result;
	}
}

static PyMethodDef PTMModuleMethods[] =
{
	{"index_structure", (PyCFunction)index_structure, METH_VARARGS | METH_KEYWORDS, "determine the structure of the atom"},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initptmmodule(void)
{
	(void) Py_InitModule("ptmmodule", PTMModuleMethods);
	import_array();

	initialize_PTM();
	//uint64_t res = run_tests();
	//if (res != 0)
	//	return error(PyExc_RuntimeError, "PTM unit tests failed");
}

#ifdef __cplusplus
}
#endif

