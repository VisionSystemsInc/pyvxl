#include "pyvnl.h"
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_quaternion.h>
#include <tuple>
#include "pyvxl_util.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>


namespace py = pybind11;

namespace pyvxl { namespace vnl {
// VNL HELPER FUNCTIONS

template <class T>
T vnl_vector_getitem(vnl_vector<T> const& v, long i)
{
  // wrap around
  if (i < 0) {
    i += v.size();
  }
  if ((i < 0) || (i >= v.size())) {
    throw py::index_error("index out of range");
  }
  return v[i];
}

template <class T>
int vnl_vector_len(vnl_vector<T> const& v)
{
  return v.size();
}

template <class T>
vnl_vector<T> vnl_matrix_getitem(vnl_matrix<T> const& m, long r)
{
  if (r < 0) {
    r += m.rows();
  }
  if ((r < 0) || (r >= m.rows())) {
    throw py::index_error("index out of range");
  }
  return m.get_row(r);
}

template<class T, unsigned R, unsigned C>
vnl_vector<T> vnl_matrix_fixed_getitem(vnl_matrix_fixed<T,R,C> const& m, long r)
{
  return vnl_matrix_getitem<T>(m.as_ref(), r);
}

template<class T, unsigned N>
T vnl_vector_fixed_getitem(vnl_vector_fixed<T,N> const& v, long i)
{
  return vnl_vector_getitem<T>(v.as_ref(), i);
}

template<class T>
T vnl_quaternion_getitem(vnl_quaternion<T> const& q, long i)
{
  // wrap around
  if (i < 0) {
    i += 4;
  }
  if ((i < 0) || (i >= 4)) {
    throw py::index_error("index out of range");
  }
  switch(i) {
    case 0: return q.x();
    case 1: return q.y();
    case 2: return q.z();
    case 3: return q.r();
  }
  return NAN;
}


template <class T>
int vnl_matrix_len(vnl_matrix<T> const& m)
{
  return m.rows();
}

template<class T, unsigned R, unsigned C>
int vnl_matrix_fixed_len(vnl_matrix_fixed<T,R,C> const& m)
{
  return R;
}

template<class T, unsigned N>
int vnl_vector_fixed_len(vnl_vector_fixed<T,N> const& v)
{
  return N;
}

template<class T>
int vnl_quaternion_len(vnl_quaternion<T> const& q)
{
  return 4;
}

template<class T>
std::tuple<int, int> vnl_matrix_shape(vnl_matrix<T> const& m)
{
  return std::make_tuple<int,int>(m.rows(), m.cols());
}

template<class T>
std::tuple<int> vnl_vector_shape(vnl_vector<T> const& v)
{
  return std::make_tuple<int>(v.size());
}

template<class T, unsigned R, unsigned C>
std::tuple<int, int> vnl_matrix_fixed_shape(vnl_matrix_fixed<T,R,C> const& m)
{
  return std::make_tuple<int,int>(R,C);
}

template<class T, unsigned N>
std::tuple<int> vnl_vector_fixed_shape(vnl_vector_fixed<T,N> const& v)
{
  return std::make_tuple<int>(N);
}

template<class T>
std::tuple<int> vnl_quaternion_shape(vnl_quaternion<T> const& q)
{
  return std::make_tuple<int>(4);
}

template<class T, unsigned NR, unsigned NC>
vnl_matrix_fixed<T,NR,NC>* matrix_fixed_from_buffer(py::array_t<T> b)
{
  py::buffer_info info = b.request();
  if (info.format != py::format_descriptor<T>::format()) {
    throw std::runtime_error("Incompatible scalar type");
  }
  if (info.ndim != 2) {
    throw std::runtime_error("Expecting a 2-dimensional matrix");
  }
  const size_t num_rows = info.shape[0];
  if (num_rows != NR) {
    std::stringstream errstr;
    errstr << "Input must have " << NR << " rows.";
    throw std::runtime_error(errstr.str().c_str());
  }
  const size_t num_cols = info.shape[1];
  if (num_cols != NC) {
    std::stringstream errstr;
    errstr << "Input must have " << NC << " cols.";
    throw std::runtime_error(errstr.str().c_str());
  }

  vnl_matrix_fixed<T,NR,NC> *mat = new vnl_matrix_fixed<T,NR,NC>();
  int row_stride = info.strides[0]/sizeof(T);
  int col_stride = info.strides[1]/sizeof(T);
  for (size_t r=0; r<NR; ++r) {
    for (size_t c=0; c<NC; ++c) {
      (*mat)(r,c) = *(static_cast<T*>(info.ptr) + r*row_stride + c*col_stride);
    }
  }
  return mat;
}

template<class T, unsigned N>
vnl_vector_fixed<T,N>* vector_fixed_from_buffer(py::array_t<T> b)
{
  py::buffer_info info = b.request();
  if (info.format != py::format_descriptor<T>::format()) {
    throw std::runtime_error("Incompatible scalar type");
  }
  if (info.ndim != 1) {
    throw std::runtime_error("Expecting a 1-dimensional vector");
  }
  const size_t num_elements = info.shape[0];
  if (num_elements != N) {
    std::stringstream errstr;
    errstr << "Input must have " << N << " elements.";
    throw std::runtime_error(errstr.str().c_str());
  }
  vnl_vector_fixed<T,N> *vec = new vnl_vector_fixed<T,N>();
  size_t stride = info.strides[0]/sizeof(T);
  for (size_t n=0; n<N; ++n) {
    (*vec)(n) = *(static_cast<T*>(info.ptr) + n*stride);
  }

  return vec;
}

template<class T>
vnl_quaternion<T>* quaternion_from_buffer(py::array_t<T> b)
{
  vnl_vector_fixed<T,4> *qvec = vector_fixed_from_buffer<T,4>(b);
  vnl_quaternion<T> *q = new vnl_quaternion<T>(*qvec);
  delete qvec;
  return q;
}

template <class T>
vnl_matrix<T>* matrix_from_buffer(py::array_t<T> b)
{
  py::buffer_info info = b.request();
  if (info.format != py::format_descriptor<T>::format()) {
    throw std::runtime_error("Incompatible scalar type");
  }
  if (info.ndim != 2) {
    throw std::runtime_error("Expecting a 2-dimensional matrix");
  }
  const size_t num_rows = info.shape[0];
  const size_t num_cols = info.shape[1];
  vnl_matrix<T> *mat = new vnl_matrix<T>(num_rows, num_cols);
  int row_stride = info.strides[0]/sizeof(T);
  int col_stride = info.strides[1]/sizeof(T);
  for (size_t r=0; r<num_rows; ++r) {
    for (size_t c=0; c<num_cols; ++c) {
      (*mat)(r,c) = *(static_cast<T*>(info.ptr) + r*row_stride + c*col_stride);
    }
  }

  return mat;
}

template <class T>
vnl_vector<T>* vector_from_buffer(py::array_t<T> b)
{
  py::buffer_info info = b.request();
  if (info.format != py::format_descriptor<T>::format()) {
    throw std::runtime_error("Incompatible scalar type");
  }
  if (info.ndim != 1) {
    throw std::runtime_error("Expecting a 1-dimensional vector");
  }
  const size_t num_elements = info.shape[0];
  // in-place constructor
  vnl_vector<T>* vec = new vnl_vector<T>(num_elements);
  int stride = info.strides[0]/sizeof(T);
  for (size_t n=0; n<num_elements; ++n) {
    (*vec)(n) = *(static_cast<T*>(info.ptr) + n*stride);
  }
  return vec;
}

template<class T>
py::buffer_info get_matrix_buffer(vnl_matrix<T> &m)
{
  return py::buffer_info(m.data_block(), sizeof(T),
                         py::format_descriptor<T>::format(),
                         2, {m.rows(), m.cols()},
                         {sizeof(T)*m.cols(), sizeof(T)});
}

template<class T, unsigned NR, unsigned NC>
py::buffer_info get_matrix_fixed_buffer(vnl_matrix_fixed<T,NR,NC> &m)
{
  return py::buffer_info(m.data_block(), sizeof(T),
                         py::format_descriptor<T>::format(),
                         2, {NR, NC},
                         {sizeof(T)*NC, sizeof(T)});
}

template<class T>
py::buffer_info get_vector_buffer(vnl_vector<T> &v)
{
  return py::buffer_info(v.data_block(), sizeof(T),
                         py::format_descriptor<T>::format(),
                         1, {v.size()},
                         {sizeof(T)});
}

template<class T, unsigned N>
py::buffer_info get_vector_fixed_buffer(vnl_vector_fixed<T,N> &v)
{
  return py::buffer_info(v.data_block(), sizeof(T),
                         py::format_descriptor<T>::format(),
                         1, {N},
                         {sizeof(T)});
}

template<class T>
void wrap_vnl_matrix(py::module &m, std::string const& class_name)
{
  py::class_<vnl_matrix<T> >(m, class_name.c_str(), py::buffer_protocol())
    .def(py::init<unsigned int,unsigned int>())
    .def(py::init(&matrix_from_buffer<T>))
    .def("get", &vnl_matrix<T>::get)
    .def_property_readonly("shape", &vnl_matrix_shape<T>)
    .def("__str__", stream2str<vnl_matrix<T> >)
    .def("__getitem__", vnl_matrix_getitem<T>)
    .def("__len__", vnl_matrix_len<T>)
    .def(py::self + vnl_matrix<T>())
    .def(py::self * vnl_vector<T>())
    .def_buffer(get_matrix_buffer<T>);

  py::implicitly_convertible<py::array_t<T>, vnl_matrix<T> >();

}

template<class T>
void wrap_vnl_vector(py::module &m, std::string const& class_name)
{
  py::class_<vnl_vector<T> >(m, class_name.c_str(), py::buffer_protocol())
    .def(py::init<size_t>())
    .def(py::init<size_t, T>())
    .def(py::init(&vector_from_buffer<T>))
    .def("get", &vnl_vector<T>::get)
    .def("size", &vnl_vector<T>::size)
    .def("__len__", vnl_vector_len<T>)
    .def("__str__", stream2str<vnl_vector<T> >)
    .def("__getitem__", vnl_vector_getitem<T>)
    .def(py::self + py::self)
    .def_buffer(get_vector_buffer<T>);

  py::implicitly_convertible<py::array_t<T>, vnl_vector<T> >();
}

template<class T, unsigned N>
void wrap_vnl_vector_fixed(py::module &m, std::string const& class_name)
{
  py::class_<vnl_vector_fixed<T,N> >(m, class_name.c_str(), py::buffer_protocol())
    .def(py::init(&vector_fixed_from_buffer<T,N>))
    .def("get", &vnl_vector_fixed<T,N>::get)
    .def_property_readonly("shape", &vnl_vector_fixed_shape<T,N>)
    .def("__str__", stream2str<vnl_vector_fixed<T,N> >)
    .def("__len__", vnl_vector_fixed_len<T,N>)
    .def("__getitem__", vnl_vector_fixed_getitem<T,N>)
    .def(py::self + vnl_vector<T>())
    .def(py::self + py::self)
    .def_buffer(get_vector_fixed_buffer<T,N>);

  py::implicitly_convertible<py::array_t<T>, vnl_vector_fixed<T,N> >();
}

template<class T, unsigned NR, unsigned NC>
void wrap_vnl_matrix_fixed(py::module &m, std::string const& class_name)
{
  py::class_<vnl_matrix_fixed<T,NR,NC> >(m, class_name.c_str(), py::buffer_protocol())
    .def(py::init(&matrix_fixed_from_buffer<T,NR,NC>))
    .def("get", &vnl_matrix_fixed<T,NR,NC>::get)
    .def_property_readonly("shape", &vnl_matrix_fixed_shape<T,NR,NC>)
    .def("__str__", stream2str<vnl_matrix_fixed<T,NR,NC> >)
    .def("__getitem__", vnl_matrix_fixed_getitem<T,NR,NC>)
    .def("__len__", vnl_matrix_fixed_len<T,NR,NC>)
    .def(py::self + py::self)
    .def(py::self * vnl_vector<T>())
    .def(py::self * vnl_vector_fixed<T,NC>())
    .def_buffer(get_matrix_fixed_buffer<T,NR,NC>);

  py::implicitly_convertible<py::array_t<T>, vnl_matrix_fixed<T,NR,NC> >();
}

template<class T>
void wrap_vnl_quaternion(py::module &m, std::string const& class_name)
{
  py::class_<vnl_quaternion<T> >(m, class_name.c_str())
    .def(py::init(&quaternion_from_buffer<T>))
    .def_property_readonly("shape", &vnl_quaternion_shape<T>)
    .def("__str__", stream2str<vnl_quaternion<T> >)
    .def("__getitem__", vnl_quaternion_getitem<T>)
    .def("__len__", [](vnl_quaternion<T> const& q){return (size_t)4;});
}


void wrap_vnl(py::module &m)
{
  wrap_vnl_vector<double>(m, "vector");
  wrap_vnl_matrix<double>(m, "matrix");
  wrap_vnl_matrix_fixed<double,3,3>(m, "matrix_fixed_3x3");
  wrap_vnl_matrix_fixed<double,3,4>(m, "matrix_fixed_3x4");
  wrap_vnl_matrix_fixed<double,4,20>(m, "matrix_fixed_4x20");
  wrap_vnl_vector_fixed<double,3>(m, "vector_fixed_3");
  wrap_vnl_vector_fixed<double,4>(m, "vector_fixed_4");
  wrap_vnl_quaternion<double>(m, "quaternion");
}
}}
