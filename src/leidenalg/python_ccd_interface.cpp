#include "python_ccd_interface.h"
#include <iostream>

#ifdef __cplusplus
extern "C"
{
#endif

PyObject* _calcCCD(PyObject* self, PyObject* args)
{
    // Parse the arguments to get the two NumPy arrays
    PyObject *py_refmat, *py_emat;
    int refRows, refCols, ematRows, ematCols;


     if (!PyArg_ParseTuple(args, "OiiOii", &py_refmat, &refRows, &refCols, 
                        &py_emat, &ematRows, &ematCols)) {
        return NULL;
    }

    // Check if the arrays have the desired format and flags
    if (!PyArray_Check(py_refmat) || !PyArray_Check(py_emat)) {
        PyErr_SetString(PyExc_TypeError, "Expected NumPy arrays");
        return NULL;
    }

    if (PyArray_TYPE(py_emat) != NPY_DOUBLE || PyArray_TYPE(py_refmat) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "Expected NumPy arrays of type 'double'");
        return NULL;
    }
    // Check if the dimensions match the provided arguments
    if (PyArray_DIM(py_refmat, 0) != refRows || PyArray_DIM(py_refmat, 1) != refCols){
        PyErr_SetString(PyExc_ValueError, "Mismatched dimensions for ref array");
        return NULL;
    } 
    if(PyArray_DIM(py_emat, 0) != ematRows || PyArray_DIM(py_emat, 1) != ematCols) {
        PyErr_SetString(PyExc_ValueError, "Mismatched dimensions for emat array");
        return NULL;
    }

    PyArrayObject* nparray_emat = (PyArrayObject*)PyArray_FROM_OTF(py_emat, NPY_DOUBLE, NPY_ARRAY_INOUT_ARRAY);
    PyArrayObject* nparray_ref = (PyArrayObject*)PyArray_FROM_OTF(py_refmat, NPY_DOUBLE, NPY_ARRAY_INOUT_ARRAY);

    // Access the data pointers
    // double* refmat_data = (double*)PyArray_DATA(py_refmat);
    // double* emat_data = (double*)PyArray_DATA(py_emat);
    double* refmat_data = (double*)PyArray_DATA(nparray_ref);
    double* emat_data = (double*)PyArray_DATA(nparray_emat);

   if (!nparray_emat || !nparray_ref) {
        PyErr_SetString(PyExc_TypeError, "Invalid NumPy array.");
        return NULL;
    }

    // Use PyArray_FROM_OTF directly for initialization
    std::vector<double> emat(emat_data, emat_data + ematRows * ematCols);
    std::vector<double> refMat(refmat_data, refmat_data + refRows * refCols);

    // Call your calculation function
    double result = ccd_utils::calcCCDsimple(refMat, refRows, emat, ematRows, ematCols, false);

    // Return the result
    return Py_BuildValue("d", result);

}

PyObject* _calcCCS(PyObject* self, PyObject* args)
{
    // Parse the arguments to get the two NumPy arrays
    PyObject *py_refmat, *py_emat, *py_subject_groups;
    int refRows, refCols, ematRows, ematCols;


     if (!PyArg_ParseTuple(args, "OiiOiiO", &py_refmat, &refRows, &refCols, 
                        &py_emat, &ematRows, &ematCols, &py_subject_groups)) {
        return NULL;
    }

    // Check if the arrays have the desired format and flags
    if (!PyArray_Check(py_refmat) || !PyArray_Check(py_emat) || !PyArray_Check(py_subject_groups)) {
        PyErr_SetString(PyExc_TypeError, "Expected NumPy arrays");
        return NULL;
    }

    if (PyArray_TYPE(py_subject_groups) != NPY_INT ) {
        PyErr_SetString(PyExc_TypeError, "Expected subject_info array of type 'int32'");
        return NULL;
    }
    if (PyArray_TYPE(py_emat) != NPY_DOUBLE || PyArray_TYPE(py_refmat) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "Expected NumPy arrays of type 'double'");
        return NULL;
    }

    // Check if the dimensions match the provided arguments
    if (PyArray_DIM(py_refmat, 0) != refRows || PyArray_DIM(py_refmat, 1) != refCols){
        PyErr_SetString(PyExc_ValueError, "Mismatched dimensions for ref array");
        return NULL;
    } 
    if(PyArray_DIM(py_emat, 0) != ematRows || PyArray_DIM(py_emat, 1) != ematCols) {
        PyErr_SetString(PyExc_ValueError, "Mismatched dimensions for emat array");
        return NULL;
    }

    PyArrayObject* nparray_emat = (PyArrayObject*)PyArray_FROM_OTF(py_emat, NPY_DOUBLE, NPY_ARRAY_INOUT_ARRAY);
    PyArrayObject* nparray_ref = (PyArrayObject*)PyArray_FROM_OTF(py_refmat, NPY_DOUBLE, NPY_ARRAY_INOUT_ARRAY);
    PyArrayObject* nparray_subject = (PyArrayObject*)PyArray_FROM_OTF(py_subject_groups, NPY_INT, NPY_ARRAY_INOUT_ARRAY);

    // Access the data pointers
    // double* refmat_data = (double*)PyArray_DATA(py_refmat);
    // double* emat_data = (double*)PyArray_DATA(py_emat);
    double* refmat_data = (double*)PyArray_DATA(nparray_ref);
    double* emat_data = (double*)PyArray_DATA(nparray_emat);
    int* subject_data = (int*)PyArray_DATA(nparray_subject);

   if (!nparray_emat || !nparray_ref || !nparray_subject) {
        PyErr_SetString(PyExc_TypeError, "Invalid NumPy array.");
        return NULL;
    }

    // Use PyArray_FROM_OTF directly for initialization
    std::vector<double> emat(emat_data, emat_data + ematRows * ematCols);
    std::vector<double> refMat(refmat_data, refmat_data + refRows * refCols);
    std::vector<int> subject_info(subject_data, subject_data + 1 * ematCols); 
    // Call your calculation function
    std::vector<double>emat_group_sumd;
    auto resultPair = ccd_utils::sumColumnsByGroup(emat, ematRows, ematCols, subject_info, emat_group_sumd, 1);
    int num_groups = resultPair.first;
    double result = ccd_utils::calcCCS(refMat, refRows, emat_group_sumd, ematRows, num_groups);
                    
    // Return the result
    return Py_BuildValue("d", result);

}


#ifdef __cplusplus
}
#endif