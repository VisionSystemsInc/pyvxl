#include "pyvbl.h"

// standard library
#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// pybind11
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

// vxl classes
#include <vbl/vbl_array_1d.h>
#include <vbl/vbl_array_2d.h>
#include <vbl/vbl_array_3d.h>

// io classes for py::pickle
#include <vbl/io/vbl_io_array_1d.h>
#include <vbl/io/vbl_io_array_2d.h>
#include <vbl/io/vbl_io_array_3d.h>

// utilities
#include "../pyvxl_util.h"


namespace py = pybind11;

namespace pyvxl { namespace vbl {

template<typename T>
void wrap_vbl_array_1d(py::module &m, std::string const& class_name)
{
  py::class_<vbl_array_1d<T> > (m, class_name.c_str(), py::buffer_protocol())
    .def(py::init<>())
    .def(py::init([](py::buffer b) {
      py::buffer_info info = b.request();

      if (info.format != py::format_descriptor<T>::format()) {
        throw std::runtime_error("Incompatible scalar type");
      }
      if (info.ndim != 1) {
        throw std::runtime_error("Expecting a 1-dimensional vector");
      }

      const T* data_ptr = static_cast<T*>(info.ptr);
      const size_t num_elements = info.shape[0];
      const size_t stride = info.strides[0] / sizeof(T);

      // construct
      const T* begin_ptr = static_cast<const T*>(data_ptr);
      const T* end_ptr = static_cast<const T*>(data_ptr + (num_elements * stride));
      return vbl_array_1d<T>(begin_ptr, end_ptr);
    }))
    .def_buffer([](vbl_array_1d<T> &a) -> py::buffer_info {
      return py::buffer_info(
        a.begin(),  /* Pointer to buffer */
        sizeof(T),  /* Size of one scalar */
        py::format_descriptor<T>::format(),  /* Python struct-style format descriptor */
        1,  /* Number of dimensions */
        { a.size() },  /* Buffer dimensions */
        { sizeof(T) } /* Strides (in bytes) for each index */
      );
    })
    .def("__len__", &vbl_array_1d<T>::size)
    .def("__getitem__", [](const vbl_array_1d<T> &a, size_t i) {
      // wrap around
      if (i < 0) {
        i += a.size();
      }
      if (i >= 0 && i < a.size())
        return a[i];
      else
        throw py::index_error();
    })
    .def("__repr__", streamToString<vbl_array_1d<T> >)
    .def(py::pickle(&vslPickleGetState<vbl_array_1d<T> >,
                    &vslPickleSetState<vbl_array_1d<T> >))
    .def(py::self == py::self)
    ;
}


template<typename T>
void wrap_vbl_array_2d(py::module &m, std::string const& class_name)
{
  py::class_<vbl_array_2d<T> > (m, class_name.c_str(), py::buffer_protocol())
    .def(py::init<>())
    .def(py::init([](py::buffer b) {
      py::buffer_info info = b.request();

      if (info.format != py::format_descriptor<T>::format()) {
        throw std::runtime_error("Incompatible scalar type");
      }
      if (info.ndim != 2) {
        throw std::runtime_error("Expecting a 2-dimensional vector");
      }

      const T* data_ptr = static_cast<T*>(info.ptr);
      const size_t num_rows = info.shape[0];
      const size_t num_cols = info.shape[1];
      const size_t row_stride = info.strides[0] / sizeof(T);
      const size_t col_stride = info.strides[1] / sizeof(T);

      // construct - copies data because I don't see a way to set the data
      // pointer directly in vbl_array_2d
      vbl_array_2d<T> arr(info.shape[0], info.shape[1]);
      for (size_t i=0; i < num_rows; i++) {
        for (size_t j=0; j < num_cols; j++) {
          arr.put(i, j, *(data_ptr + (i * row_stride) + (j * col_stride)));
        }
      }

      return arr;
    }))
    .def_buffer([](vbl_array_2d<T> &a) -> py::buffer_info {
      return py::buffer_info(
        a.begin(),  /* Pointer to buffer */
        sizeof(T),  /* Size of one scalar */
        py::format_descriptor<T>::format(),  /* Python struct-style format descriptor */
        2,  /* Number of dimensions */
        { a.rows(), a.cols() },  /* Buffer dimensions */
        { sizeof(T) * a.cols(), sizeof(T) } /* Strides (in bytes) for each index */
      );
    })
    .def("__len__", &vbl_array_2d<T>::size)
    .def("__getitem__", [](const vbl_array_2d<T> &a, std::tuple<size_t, size_t> index) {
      size_t i = std::get<0>(index);
      size_t j = std::get<1>(index);

      // wrap around
      if (i < 0) {
        i += a.rows();
      }
      if (j < 0) {
        j += a.cols();
      }

      if (i >= 0 && i < a.rows() &&
          j >= 0 && j < a.cols())
        return a.get(i, j);
      else
        throw py::index_error();
    })
    .def("__repr__", streamToString<vbl_array_2d<T> >)
    .def(py::pickle(&vslPickleGetState<vbl_array_2d<T> >,
                    &vslPickleSetState<vbl_array_2d<T> >))
    .def(py::self == py::self)
    ;
}


void wrap_vbl(py::module &m)
{
  wrap_vbl_array_1d<double>(m, "array_1d");
  wrap_vbl_array_1d<float>(m, "array_1d_float");

  wrap_vbl_array_2d<double>(m, "array_2d");
  wrap_vbl_array_2d<float>(m, "array_2d_float");
}
}}

PYBIND11_MODULE(_vbl, m)
{
  m.doc() =  "Python bindings for the VBL computer vision libraries";

  pyvxl::vbl::wrap_vbl(m);
}
