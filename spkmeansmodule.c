#include <Python.h>
#include "spkmeans.h"

double **mat_from_python_to_C(PyObject *_PyData, int lines, int columns)
{
    PyObject *PyRow;
    int i, j;
    double coordinate;
    double **mat = (double **)malloc(lines*sizeof(double *));
    double *mat_inner = (double *)malloc(lines*columns*sizeof(double));
    assert(mat != NULL);
    assert(mat_inner != NULL);

    for (i = 0; i < lines; ++i)
    {
        mat[i] = mat_inner + i * columns;
    }

    for (i = 0; i < lines; i++)
    {
        PyRow = PyList_GetItem(_PyData, i);
        if (!PyList_Check(PyRow))
        {
            return NULL;
        }
        for (j = 0; j < columns; j++)
        {
            if (!PyFloat_Check(PyList_GetItem(PyRow, j)))
            {
                return NULL;
            }
            coordinate = PyFloat_AsDouble(PyList_GetItem(PyRow, j));
            mat[i][j] = coordinate;
        }
    }
    return mat;
}

static PyObject *construct_python_matrix(double **C_mat, int lines, int columns)
{
    PyObject *_mat, *matRow;
    int i, j;
    _mat = PyList_New(lines);
    for (i = 0; i < lines; i++)
    {
        matRow = PyList_New(columns);
        for (j = 0; j < columns; j++)
        {
            PyList_SET_ITEM(matRow, j, PyFloat_FromDouble(C_mat[i][j]));
        }
        PyList_SET_ITEM(_mat, i, matRow);
    }
    return _mat;
}

static PyObject *goal_switch(PyObject *self, PyObject *args)
{
    PyObject *_PyData, *_mat, *matRow;
    Py_ssize_t Py_goal, Py_n, Py_k, Py_d;
    int i, j, goal, n, k, d;
    double **points, **wam, **ddg, **lnorm, **T;
    EIGEN *eigens;
    if (!PyArg_ParseTuple(args, "Oiiii", &_PyData, &Py_goal, &Py_n, &Py_k, &Py_d))
    {
        return NULL;
    }

    if (!PyList_Check(_PyData) && !PyLong_Check(Py_goal) && !PyLong_Check(Py_n) && !PyLong_Check(Py_k) && !PyLong_Check(Py_d))
    {
        return NULL;
    }
    
    goal = (int)Py_goal;
    n = (int)Py_n;
    k = (int)Py_k;
    d = (int)Py_d;

    points = mat_from_python_to_C(_PyData, n, d);

    if (goal == 1)
    { // spk
        if (k >= n)
        {
            printf("Invalid Input!\n");
            exit(0);
        }
        wam = gen_wam(points, n, d);
        ddg = gen_ddg(wam, n);
        lnorm = gen_lnorm(wam, ddg, n);
        eigens = jacobi(lnorm, n);
        qsort(eigens, n, sizeof(EIGEN), my_comparator);
        if (k == 0)
        {
            k = determine_k(eigens, n);
        }
        d = k;
        T = gen_T(eigens, n, k);
        _mat = construct_python_matrix(T, n, k);

        free(wam[0]);
        free(wam);
        free(ddg[0]);
        free(ddg);
        free(lnorm[0]);
        free(lnorm);
        for (i = 0; i < n; i++)
        {
            free(eigens[i].vector);
        }
        free(eigens);
        free(points[0]);
        free(points);
        free(T[0]);
        free(T);
        
        return _mat;
    }
    else if (goal == 2)
    { // wam
        wam = gen_wam(points, n, d);
        _mat = construct_python_matrix(wam, n, n);

        free(wam[0]);
        free(wam);

        return _mat;
    }
    else if (goal == 3)
    { // ddg
        wam = gen_wam(points, n, d);
        ddg = gen_ddg(wam, n);
        _mat = construct_python_matrix(ddg, n, n);

        free(wam[0]);
        free(wam);
        free(ddg[0]);
        free(ddg);

        return _mat;
    }
    else if (goal == 4)
    { // lnorm
        wam = gen_wam(points, n, d);
        ddg = gen_ddg(wam, n);
        lnorm = gen_lnorm(wam, ddg, n);
        _mat = construct_python_matrix(lnorm, n, n);

        free(wam[0]);
        free(wam);
        free(ddg[0]);
        free(ddg);
        free(lnorm[0]);
        free(lnorm);

        return _mat;
    }
    else
    { //goal == 5, jacobi
        wam = gen_wam(points, n, d);
        ddg = gen_ddg(wam, n);
        lnorm = gen_lnorm(wam, ddg, n);
        eigens = jacobi(lnorm, n);

        _mat = PyTuple_New(n + 1);
        matRow = PyTuple_New(n);
        for (i = 0; i < n; i++)
        {
            PyTuple_SET_ITEM(matRow, i, PyFloat_FromDouble(eigens[i].value));
        }
        PyTuple_SET_ITEM(_mat, 0, matRow);
        for (i = 0; i < n; i++)
        {
            matRow = PyTuple_New(n);
            for (j = 0; j < n; j++)
            {
                PyTuple_SET_ITEM(matRow, j, PyFloat_FromDouble(eigens[i].vector[j]));
            }
            PyTuple_SET_ITEM(_mat, i+1, matRow);
        }

        free(wam[0]);
        free(wam);
        free(ddg[0]);
        free(ddg);
        free(lnorm[0]);
        free(lnorm);
        for (i = 0; i < n; i++)
        {
            free(eigens[i].vector);
        }
        free(eigens);

        return _mat;
    }
}

static PyObject *fit(PyObject *self, PyObject *args)
{
    PyObject *_PyData, *_PyCenters, *_ret, *retCenter;
    Py_ssize_t Py_k, Py_max_iter, Py_d, Py_n;
    int i, j, k, max_iter, d, n, changed;
    double **data, **centroids, **sums;
    int *sizes;

    if (!PyArg_ParseTuple(args, "OOiiii", &_PyData, &_PyCenters, &Py_n, &Py_k, &Py_d, &Py_max_iter))
    {
        return NULL;
    }

    if (!PyList_Check(_PyData) && !PyList_Check(_PyCenters) && !PyLong_Check(Py_k) && !PyLong_Check(Py_max_iter) && !PyLong_Check(Py_d) && !PyLong_Check(Py_n))
    {
        return NULL;
    }

    n = (int)Py_n;
    k = (int)Py_k;
    d = (int)Py_d;
    max_iter = (int)Py_max_iter;

    data = mat_from_python_to_C(_PyData, n, d);
    centroids = mat_from_python_to_C(_PyCenters, k, d);
    sizes = initialize_sizes(k);
    sums = initialize_sums(k, d);

    changed = 1;
    while (changed == 1 && max_iter > 0)
    {
        assign_vec_to_cluster(data, sizes, sums, centroids, n, k, d);
        changed = update_clusters(sizes, sums, centroids, k, d);
        max_iter -= 1;
    }

    free(data[0]);
    free(data);

    free(sums[0]);
    free(sums);

    free(sizes);

    _ret = PyTuple_New(k);
    for (i = 0; i < k; i++)
    {
        retCenter = PyTuple_New(d);
        for (j = 0; j < d; j++)
        {
            PyTuple_SET_ITEM(retCenter, j, PyFloat_FromDouble(centroids[i][j]));
        }
        PyTuple_SET_ITEM(_ret, i, retCenter);
    }

    free(centroids[0]);
    free(centroids);

    return _ret;
}

static PyMethodDef capiMethods[] = {
    {"goal_switch", (PyCFunction)goal_switch,METH_VARARGS, PyDoc_STR("output matrix according goal")},
    {"fit", (PyCFunction)fit, METH_VARARGS, PyDoc_STR("perform the last step in spk")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule",
    NULL,
    -1,
    capiMethods};

PyMODINIT_FUNC
PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}