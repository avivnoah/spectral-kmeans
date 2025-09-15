# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include "symnmf.h"

static PyObject* sym(PyObject *self, PyObject *args){
    PyObject *elements_seq = NULL;
    PyObject *element = NULL;
    PyObject *item = NULL;
    double element_coord;
    double **elements;
    double **returned_mat;
    int i, j;
    int num_of_elements;
    int d;

    int argc = PyTuple_Size(args);
    if (argc == 3) {
        if(!PyArg_ParseTuple(args, "Oii", &elements_seq,  &num_of_elements, &d)) {
            return NULL;
        }
    } else if (argc == 1) {
        if(!PyArg_ParseTuple(args, "O", &elements_seq)) {
            return NULL;
        }
        num_of_elements = (int)PySequence_Length(elements_seq);
        if (num_of_elements <= 0) {
            PyErr_SetString(PyExc_ValueError, "empty input");
            return NULL;
        }
        element = PySequence_GetItem(elements_seq, 0); /* new ref */
        if (!element) return NULL;
        d = (int)PySequence_Length(element);
        Py_DECREF(element);
    } else {
        PyErr_SetString(PyExc_TypeError, "sym expects 1 or 3 arguments");
        return NULL;
    }

    /*memory allocation for all the points in the file*/
    elements = matrix_allocation(num_of_elements, d);

    /*reading the data points from python and passing into c matrix*/
    for(i=0; i<num_of_elements; i++){
        element = PySequence_GetItem(elements_seq, i); /* new ref */
        if (!element) {
            free_matrix(elements);
            PyErr_SetString(PyExc_IndexError, "failed to get element row");
            return NULL;
        }
        /* validate inner length */
        if ((int)PySequence_Length(element) != d) {
            Py_DECREF(element);
            free_matrix(elements);
            PyErr_SetString(PyExc_ValueError, "inconsistent inner sequence length");
            return NULL;
        }
        for(j=0; j<d; j++){
            item = PySequence_GetItem(element, j); /* new ref */
            if (!item) {
                Py_DECREF(element);
                free_matrix(elements);
                PyErr_SetString(PyExc_IndexError, "failed to get coordinate");
                return NULL;
            }
            element_coord = PyFloat_AsDouble(item);
            Py_DECREF(item);
            if (PyErr_Occurred()) {
                Py_DECREF(element);
                free_matrix(elements);
                return NULL;
            }
            elements[i][j] = element_coord;
        }
        Py_DECREF(element);
    }

    returned_mat = sym_c(elements, num_of_elements, d);

    free_matrix(elements);

    PyObject* matrix;
    PyObject* vector;
    PyObject* python_float;
    matrix = PyList_New(num_of_elements);
    for (int ii = 0; ii < num_of_elements; ii++)
    {   
        vector = PyList_New(num_of_elements);
        for(j=0; j<num_of_elements; j++){
            python_float = PyFloat_FromDouble(returned_mat[ii][j]);
            PyList_SetItem(vector, j, python_float); 
        }
        PyList_SetItem(matrix, ii, vector);
    }

    free_matrix(returned_mat);

    return Py_BuildValue("O", matrix); 
}

static PyObject* ddg(PyObject *self, PyObject *args){
    PyObject *elements_seq = NULL;
    PyObject *element = NULL;
    PyObject *item = NULL;
    double element_coord;
    double **elements;
    double **returned_mat;
    int i, j;
    int num_of_elements;
    int d;

    int argc = PyTuple_Size(args);
    if (argc == 3) {
        if(!PyArg_ParseTuple(args, "Oii", &elements_seq,  &num_of_elements, &d)) {
            return NULL;
        }
    } else if (argc == 1) {
        if(!PyArg_ParseTuple(args, "O", &elements_seq)) {
            return NULL;
        }
        num_of_elements = (int)PySequence_Length(elements_seq);
        if (num_of_elements <= 0) {
            PyErr_SetString(PyExc_ValueError, "empty input");
            return NULL;
        }
        element = PySequence_GetItem(elements_seq, 0);
        if (!element) return NULL;
        d = (int)PySequence_Length(element);
        Py_DECREF(element);
    } else {
        PyErr_SetString(PyExc_TypeError, "ddg expects 1 or 3 arguments");
        return NULL;
    }

    /*memory allocation for all the points in the file*/
    elements = matrix_allocation(num_of_elements, d);

    /*reading the data points from python and passing into c matrix*/
    for(i=0; i<num_of_elements; i++){
        element = PySequence_GetItem(elements_seq, i);
        if (!element) {
            free_matrix(elements);
            PyErr_SetString(PyExc_IndexError, "failed to get element row");
            return NULL;
        }
        /* validate inner length */
        if ((int)PySequence_Length(element) != d) {
            Py_DECREF(element);
            free_matrix(elements);
            PyErr_SetString(PyExc_ValueError, "inconsistent inner sequence length");
            return NULL;
        }
        for(j=0; j<d; j++){
            item = PySequence_GetItem(element, j);
            if (!item) {
                Py_DECREF(element);
                free_matrix(elements);
                PyErr_SetString(PyExc_IndexError, "failed to get coordinate");
                return NULL;
            }
            element_coord = PyFloat_AsDouble(item);
            Py_DECREF(item);
            if (PyErr_Occurred()) {
                Py_DECREF(element);
                free_matrix(elements);
                return NULL;
            }
            elements[i][j] = element_coord;
        }
        Py_DECREF(element);
    }

    returned_mat = ddg_c(elements, num_of_elements, d);

    free_matrix(elements);

    PyObject* matrix;
    PyObject* vector;
    PyObject* python_float;
    matrix = PyList_New(num_of_elements);
    for (int ii = 0; ii < num_of_elements; ii++){   
        vector = PyList_New(num_of_elements);
        for(j=0; j<num_of_elements; j++){
            python_float = PyFloat_FromDouble(returned_mat[ii][j]);
            PyList_SetItem(vector, j, python_float); 
        }
        PyList_SetItem(matrix, ii, vector);
    }

    free_matrix(returned_mat);

    return Py_BuildValue("O", matrix); 
}

static PyObject* norm(PyObject *self, PyObject *args){
    PyObject *elements_seq = NULL;
    PyObject *element = NULL;
    PyObject *item = NULL;
    double element_coord;
    double **elements;
    double **returned_mat;
    int i, j;
    int num_of_elements;
    int d;

    int argc = PyTuple_Size(args);
    if (argc == 3) {
        if(!PyArg_ParseTuple(args, "Oii", &elements_seq,  &num_of_elements, &d)) {
            return NULL;
        }
    } else if (argc == 1) {
        if(!PyArg_ParseTuple(args, "O", &elements_seq)) {
            return NULL;
        }
        num_of_elements = (int)PySequence_Length(elements_seq);
        if (num_of_elements <= 0) {
            PyErr_SetString(PyExc_ValueError, "empty input");
            return NULL;
        }
        element = PySequence_GetItem(elements_seq, 0);
        if (!element) return NULL;
        d = (int)PySequence_Length(element);
        Py_DECREF(element);
    } else {
        PyErr_SetString(PyExc_TypeError, "norm expects 1 or 3 arguments");
        return NULL;
    }

    /*memory allocation for all the points in the file*/
    elements = matrix_allocation(num_of_elements, d);

    /*reading the data points from python and passing into c matrix*/
    for(i=0; i<num_of_elements; i++){
        element = PySequence_GetItem(elements_seq, i);
        if (!element) {
            free_matrix(elements);
            PyErr_SetString(PyExc_IndexError, "failed to get element row");
            return NULL;
        }
        /* validate inner length */
        if ((int)PySequence_Length(element) != d) {
            Py_DECREF(element);
            free_matrix(elements);
            PyErr_SetString(PyExc_ValueError, "inconsistent inner sequence length");
            return NULL;
        }
        for(j=0; j<d; j++){
            item = PySequence_GetItem(element, j);
            if (!item) {
                Py_DECREF(element);
                free_matrix(elements);
                PyErr_SetString(PyExc_IndexError, "failed to get coordinate");
                return NULL;
            }
            element_coord = PyFloat_AsDouble(item);
            Py_DECREF(item);
            if (PyErr_Occurred()) {
                Py_DECREF(element);
                free_matrix(elements);
                return NULL;
            }
            elements[i][j] = element_coord;
        }
        Py_DECREF(element);
    }

    returned_mat = norm_c(elements, num_of_elements, d);

    free_matrix(elements);

    PyObject* matrix;
    PyObject* vector;
    PyObject* python_float;
    matrix = PyList_New(num_of_elements);
    for (int ii = 0; ii < num_of_elements; ii++)
    {   
        vector = PyList_New(num_of_elements);
        for(j=0; j<num_of_elements; j++){
            python_float = PyFloat_FromDouble(returned_mat[ii][j]);
            PyList_SetItem(vector, j, python_float); 
        }
        PyList_SetItem(matrix, ii, vector);
    }

    free_matrix(returned_mat);

    return Py_BuildValue("O", matrix); 
}

static PyObject* symnmf(PyObject *self, PyObject *args){
    PyObject *H_obj = NULL;
    PyObject *W_obj = NULL;
    PyObject *H = NULL;
    PyObject *W = NULL;
    PyObject *H_item = NULL;
    PyObject *W_item = NULL;
    PyObject *python_float = NULL;
    double H_entry;
    double W_entry;
    double **H_c;
    double **W_c;
    double **returned_mat;
    int i, j;
    int num_of_elements = -1;
    int k = -1;
    int max_iter = 300;
    /*double eps = 1e-4;*/

    int argc = PyTuple_Size(args);
    if (argc == 4) {
        if(!PyArg_ParseTuple(args, "OOii", &H_obj, &W_obj, &k, &num_of_elements)) {
            return NULL;
        }
    } else if (argc == 2) {
        if(!PyArg_ParseTuple(args, "OO", &H_obj, &W_obj)) {
            return NULL;
        }
        /* infer sizes from sequences */
        num_of_elements = (int)PySequence_Length(H_obj);
        if (num_of_elements <= 0) { PyErr_SetString(PyExc_ValueError, "empty H"); return NULL; }
        H_item = PySequence_GetItem(H_obj, 0);
        if (!H_item) return NULL;
        k = (int)PySequence_Length(H_item);
        Py_DECREF(H_item);
    } else {
        PyErr_SetString(PyExc_TypeError, "symnmf expects 2 or 4 arguments");
        return NULL;
    }

    /* allocate H_c and W_c */
    H_c = matrix_allocation(num_of_elements, k);
    W_c = matrix_allocation(num_of_elements, num_of_elements);

    /*reading H from python and passing into c matrix*/
    for(i=0; i<num_of_elements; i++){
        H = PySequence_GetItem(H_obj, i);
        if (!H) {
            free_matrix(H_c);
            free_matrix(W_c);
            PyErr_SetString(PyExc_IndexError, "failed to get H row");
            return NULL;
        }
        /* validate inner length */
        if ((int)PySequence_Length(H) != k) {
            Py_DECREF(H);
            free_matrix(H_c);
            free_matrix(W_c);
            PyErr_SetString(PyExc_ValueError, "inconsistent H inner length");
            return NULL;
        }
        for(j=0; j<k; j++){
            H_item = PySequence_GetItem(H, j);
            if (!H_item) {
                Py_DECREF(H);
                free_matrix(H_c);
                free_matrix(W_c);
                PyErr_SetString(PyExc_IndexError, "failed to get H entry");
                return NULL;
            }
            H_entry = PyFloat_AsDouble(H_item);
            Py_DECREF(H_item);
            if (PyErr_Occurred()) {
                Py_DECREF(H);
                free_matrix(H_c);
                free_matrix(W_c);
                return NULL;
            }
            H_c[i][j] = H_entry;
        }
        Py_DECREF(H);
    }

    /*reading W from python and passing into c matrix*/
    for(i=0; i<num_of_elements; i++){
        W = PySequence_GetItem(W_obj, i);
        if (!W) {
            free_matrix(H_c);
            free_matrix(W_c);
            PyErr_SetString(PyExc_IndexError, "failed to get W row");
            return NULL;
        }
        /* validate inner length */
        if ((int)PySequence_Length(W) != num_of_elements) {
            Py_DECREF(W);
            free_matrix(H_c);
            free_matrix(W_c);
            PyErr_SetString(PyExc_ValueError, "inconsistent W inner length");
            return NULL;
        }
        for(j=0; j<num_of_elements; j++){
            W_item = PySequence_GetItem(W, j);
            if (!W_item) {
                Py_DECREF(W);
                free_matrix(H_c);
                free_matrix(W_c);
                PyErr_SetString(PyExc_IndexError, "failed to get W entry");
                return NULL;
            }
            W_entry = PyFloat_AsDouble(W_item);
            Py_DECREF(W_item);
            if (PyErr_Occurred()) {
                Py_DECREF(W);
                free_matrix(H_c);
                free_matrix(W_c);
                return NULL;
            }
            W_c[i][j] = W_entry;
        }
        Py_DECREF(W);
    }

    returned_mat = symnmf_c(H_c, W_c, k, num_of_elements, max_iter);

    free_matrix(H_c);
    free_matrix(W_c);

    PyObject* matrix;
    PyObject* vector;
    matrix = PyList_New(num_of_elements);
    for (int ii = 0; ii < num_of_elements; ii++)
    {   
        vector = PyList_New(k);
        for(j=0; j<k; j++){
            python_float = PyFloat_FromDouble(returned_mat[ii][j]);
            PyList_SetItem(vector, j, python_float); 
        }
        PyList_SetItem(matrix, ii, vector);
    }

    free_matrix(returned_mat);

    return Py_BuildValue("O", matrix);
}

static PyMethodDef symnmfMethods[] = {
    {"sym",                   /* the Python method name that will be used */
      (PyCFunction) sym, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("Calculates and outputs the similarity martix from given data points.\n"
                "The similarity matrix A has num_of_elements rows and num_of_elements columns\n"
                "a_ij = exp(-0.5*(Euclidean distance(x_i-x_j))^2) if i!=j, or 0 if i=j \n"
                "expected arguments: \n"
                "X- A 2D list of d dimentional data points, of type float. Denoted by x_1, x_2,...\n"
                "num_of_elements- The number of points in X. int.\n"
                "d- The number of coordinates of each point. int.")}, /*  The docstring for the function */
      {"ddg",                   /* the Python method name that will be used */
      (PyCFunction) ddg, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("Calculates and outputs the diagonal degree martix from given data points.\n"
                "The diagonal degree martix D has num_of_elements rows and num_of_elements columns\n"
                "d_ii equals to the sum of the i'th row of the similarity matrix, and 0 elsewhwre. \n"
                "expected arguments: \n"
                "X- A 2D list of d dimentional data points, of type float. Denoted by x_1, x_2,...\n"
                "num_of_elements- the number of points in X. int.\n"
                "d- number of coordinates of each point. int.")}, /*  The docstring for the function */
      {"norm",                   /* the Python method name that will be used */
      (PyCFunction) norm, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("Calculates and outputs the normalized similarity martix W from given data points.\n"
                "The normalized similarity martix W has num_of_elements rows and num_of_elements columns\n"
                "W = D^-0.5*A*D^-0.5 where D is the diagonal degree martix, and A is the similariry matrix.\n"
                "expected arguments: \n"
                "X- A 2D list of d dimentional data points, of type float. Denoted by x_1, x_2,...\n"
                "num_of_elements- the number of points in X. int.\n"
                "d- number of coordinates of each point. int.")}, /*  The docstring for the function */
      {"symnmf",                   /* the Python method name that will be used */
      (PyCFunction) symnmf, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("Performs full the symNMF algorithm and output the final H.\n"
                "The output matrix H has num_of_elements rows and k columns.\n"
                "expected arguments: \n"
                "H- A non-negative matrix of num_of_elements rows by k columns. \n"
                "H is randomly initialized with values from the interval [0; 2*sqrt(m/k)], \n"
                "where m is the average of all entries of W. Type float.\n"
                "W- The normalized similarity matrix. Created from num_of_elements data points. Type float.\n"
                "k- the number of required clusters. int.\n"
                "num_of_elements- the number of data points that were used to create W. int.\n")}, /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL */
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf_c", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    symnmfMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_symnmf_c(void){
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}