#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <math.h>
#include "kmeans.h"
#include "struct.h"

static PyObject* fit(PyObject *self, PyObject*args) {
    /*declaring variables*/
    int K, iter, d, num_points, i, j;
    double eps, val_coord;
    PyObject* centroids_lst_Py;
    PyObject* arr_points_Py;
    Point *all_points;
    PyObject* item_point;
    PyObject* coord;
    PyObject* centroid;
    Centroid *centroids_lst;

    if (!PyArg_ParseTuple(args, "OOiiiid", &arr_points_Py, &centroids_lst_Py, &K, &iter, &d, &num_points, &eps)) {
        return NULL;
    }
    all_points = (Point*)malloc(num_points * sizeof(Point));
    centroids_lst = (Centroid*)malloc(K * sizeof(Centroid));

    /*check if allocation was successful*/
    if (all_points == NULL) {
        return NULL;
    }
    if (centroids_lst == NULL) {
        return NULL;
    }
    /*allocate space for both points & centroids*/
    for (j=0; j < num_points; j++) {
        all_points[j].coordinates = (double *) malloc(d* sizeof(double));
    }

    for (i=0; i < K; i++) {
        centroids_lst[i].coordinates = (double *) malloc(d* sizeof(double));
        centroids_lst[i].numPoints = 0;
        centroids_lst[i].Points = NULL;
    }

    /*insert arr points from python's list to all points*/
    for (j=0; j < num_points; j++) {
        item_point = PyList_GetItem(arr_points_Py,j);
        for (i = 0; i < d; i++) {
            coord = PyList_GetItem(item_point, i);
            val_coord = PyFloat_AsDouble(coord);
            all_points[j].coordinates[i] = val_coord;
        }
    }

    /*insert centroids from python's list to centroids_lst*/
    for(j = 0; j < K; j++) {
        centroid = PyList_GetItem(centroids_lst_Py,j);
        for (i = 0; i < d; i++) {
            coord = PyList_GetItem(centroid, i);
            val_coord = PyFloat_AsDouble(coord);
            centroids_lst[j].coordinates[i] = val_coord;
        }
    }
    /*Run K means algorithm*/
    kmeans_algorithm(all_points, centroids_lst, K, d, iter, num_points, eps);

    /*Convert return values to Python*/
    for(i = 0; i < K; i++) {
        item_point = PyList_GetItem(centroids_lst_Py, i);
        for (j = 0; j < d; j++) {
            coord = PyFloat_FromDouble(centroids_lst[i].coordinates[j]);
            PyList_SetItem(item_point, j, coord);
        }
    }

    /*Free memory*/
    for (j = 0; j < num_points; j++) {
        free(all_points[j].coordinates);
    }
    free(all_points);

    for (j = 0; j < K; j++) {
        free(centroids_lst[j].coordinates);
        free(centroids_lst[j].Points);
    }    
    free(centroids_lst);

    Py_INCREF(centroids_lst_Py);
    return centroids_lst_Py;
}

static PyMethodDef kmeansMethods[] = {
    {"fit", fit, METH_VARARGS, "Computes K-Means clustering."},
    {NULL, NULL, 0, NULL} 
};


static struct PyModuleDef kmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanspp", 
    NULL, 
    -1, 
    kmeansMethods 
};


PyMODINIT_FUNC PyInit_mykmeanspp(void) {
    PyObject *m;
    m = PyModule_Create(&kmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}