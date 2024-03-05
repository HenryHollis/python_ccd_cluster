#ifndef PYNTERFACE_CCD_H_INCLUDED
#define PYNTERFACE_CCD_H_INCLUDED
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL my_ARRAY_API
#include <Python.h>
#include <numpy/arrayobject.h>
#include "ccd_utils.h"

#ifdef __cplusplus
extern "C"
{
#endif

PyObject* _calcCCD(PyObject* self, PyObject* args);
#ifdef __cplusplus
}
#endif

#endif // PYNTERFACE_CCD_H_INCLUDED 
