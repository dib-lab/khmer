#include "Python.h"
#include <exception>
#include <string>
#include "oxli/oxli_exception.hh"
#include "oxli_exception_convert.hh"


void oxli_raise_py_error()
{
  try {
    throw;
  }
  catch (oxli::ReadOnlyAttribute& e) {
    PyErr_SetString(PyExc_AttributeError, e.what());
  }
  catch (oxli::InvalidValue& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
  }
  catch (oxli::InvalidStream& e) {
    PyErr_SetString(PyExc_OSError, e.what());
  }
  catch (oxli::oxli_value_exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
  }
  catch (oxli::oxli_file_exception& e) {
    PyErr_SetString(PyExc_OSError, e.what());
  }
  catch (oxli::oxli_exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
  }
}
