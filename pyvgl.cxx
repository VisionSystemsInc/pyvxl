#include "pyvgl.h"
#include <vgl/vgl_vector_2d.h>
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_point_3d.h>
#include <vgl/algo/vgl_rotation_3d.h>
#include <vgl/vgl_pointset_3d.h>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_cylinder.h>
#include <vgl/vgl_sphere_3d.h>
#include <vgl/vgl_polygon.h>
#include <vgl/vgl_line_segment_2d.h>
#include <vgl/vgl_line_segment_3d.h>
#include <vgl/vgl_oriented_box_2d.h>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <ios>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

namespace pyvxl { namespace vgl {

// Helper functions

template<class T>
double getitem_3d(T const& a, long i)
{
  // wrap around
  if (i < 0) {
    i += 3;
  }
  if (i==0) {
    return a.x();
  }
  else if (i==1) {
    return a.y();
  }
  else if (i==2) {
    return a.z();
  }
  else {
    throw py::index_error("index out of range");
  }
  return 0; // to avoid compiler warning
}

template<class T>
double getitem_2d(T const& a, long i)
{
  // wrap around
  if (i < 0) {
    i += 2;
  }
  if (i==0) {
    return a.x();
  }
  else if (i==1) {
    return a.y();
  }
  else {
    throw py::index_error("index out of range");
  }
  return 0; // to avoid compiler warning
}

template <class T, class BUFF_T>
T type_from_buffer_2d(py::array_t<BUFF_T> b)
{
  py::buffer_info info = b.request();
  if (info.format != py::format_descriptor<BUFF_T>::format()) {
    throw std::runtime_error("Incompatible scalar type");
  }
  if (info.ndim != 1) {
    throw std::runtime_error("Expecting a 1-dimensional vector");
  }
  const size_t num_elements = info.shape[0];
  if (num_elements != 2) {
    throw std::runtime_error("Expecting 2-d input vector");
  }
  // in-place constructor
  const BUFF_T* data_ptr = static_cast<BUFF_T*>(info.ptr);
  const size_t stride = info.strides[0] / sizeof(BUFF_T);
  BUFF_T x = *data_ptr;
  BUFF_T y = *(data_ptr + stride);

  return T(x,y);
}

template <class T, class BUFF_T>
T type_from_buffer_3d(py::array_t<BUFF_T> b)
{
  py::buffer_info info = b.request();
  if (info.format != py::format_descriptor<BUFF_T>::format()) {
    throw std::runtime_error("Incompatible scalar type");
  }
  if (info.ndim != 1) {
    throw std::runtime_error("Expecting a 1-dimensional vector");
  }
  const size_t num_elements = info.shape[0];
  if (num_elements != 3) {
    throw std::runtime_error("Expecting 3-d input vector");
  }
  // in-place constructor
  const BUFF_T* data_ptr = static_cast<BUFF_T*>(info.ptr);
  const size_t stride = info.strides[0] / sizeof(BUFF_T);
  BUFF_T x = *data_ptr;
  BUFF_T y = *(data_ptr + stride);
  BUFF_T z = *(data_ptr + 2*stride);

  return T(x,y,z);
}

template<typename T>
typename vgl_polygon<T>::sheet_t getitem_sheet(vgl_polygon<T> const& p, long i){

  // wrap around
  if(i < 0){
    i += p.num_sheets();
  }
  
  // out of range
  if(i < 0 || i >= p.num_sheets()){
    throw py::index_error("index out of range");
  }

  // extract sheet
  return p[i];
}

template<typename T>
std::string streamToString(T const& t){

  std::ostringstream buffer;
  buffer << t;
  return buffer.str();
}

void wrap_vgl(py::module &m)
{
  py::class_<vgl_point_2d<double> > (m, "point_2d")
    .def(py::init<double,double>())
    .def(py::init(&type_from_buffer_2d<vgl_point_2d<double>, double>))
    .def("__len__", [](vgl_point_2d<double>){return (size_t)2;})
    .def("__getitem__", getitem_2d<vgl_point_2d<double> >)
    .def("__repr__", streamToString<vgl_point_2d<double> >)
    .def_property_readonly("x", (double (vgl_point_2d<double>::*)() const) &vgl_point_2d<double>::x)
    .def_property_readonly("y", (double (vgl_point_2d<double>::*)() const) &vgl_point_2d<double>::y)
    .def(py::self - py::self)
    .def(py::self == py::self);

  py::class_<vgl_vector_2d<double> > (m, "vector_2d")
    .def(py::init<double,double>())
    .def(py::init(&type_from_buffer_2d<vgl_vector_2d<double>, double>))
    .def("__len__", [](vgl_vector_2d<double>){return (size_t)2;})
    .def("__getitem__",getitem_2d<vgl_vector_2d<double> >)
    .def("__repr__", streamToString<vgl_vector_2d<double> >)
    .def_property_readonly("x", &vgl_vector_2d<double>::x)
    .def_property_readonly("y", &vgl_vector_2d<double>::y)
    .def("length", &vgl_vector_2d<double>::length)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self == py::self);

  py::class_<vgl_point_3d<double> > (m, "point_3d")
    .def(py::init<double,double,double>())
    .def(py::init(&type_from_buffer_3d<vgl_point_3d<double>, double>))
    .def("__len__", [](vgl_point_3d<double>){return (size_t)3;})
    .def("__getitem__", getitem_3d<vgl_point_3d<double> >)
    .def("__repr__", streamToString<vgl_point_3d<double> >)
    .def_property_readonly("x", (double (vgl_point_3d<double>::*)() const) &vgl_point_3d<double>::x)
    .def_property_readonly("y", (double (vgl_point_3d<double>::*)() const) &vgl_point_3d<double>::y)
    .def_property_readonly("z", (double (vgl_point_3d<double>::*)() const) &vgl_point_3d<double>::z)
    .def(py::self - py::self)
    .def(py::self == py::self);

  py::class_<vgl_vector_3d<double> > (m, "vector_3d")
    .def(py::init<double,double,double>())
    .def(py::init(&type_from_buffer_3d<vgl_vector_3d<double>, double>))
    .def("__len__", [](vgl_vector_3d<double>){return (size_t)3;})
    .def("__getitem__", getitem_3d<vgl_vector_3d<double> >)
    .def("__repr__", streamToString<vgl_vector_3d<double> >)
    .def_property_readonly("x", &vgl_vector_3d<double>::x)
    .def_property_readonly("y", &vgl_vector_3d<double>::y)
    .def_property_readonly("z", &vgl_vector_3d<double>::z)
    .def("length", &vgl_vector_3d<double>::length)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self == py::self);

  py::class_ <vgl_rotation_3d<double> > (m, "rotation_3d")
    .def(py::init<vnl_vector_fixed<double,4> >())
    .def(py::init<vnl_matrix_fixed<double,3,3> >())
    .def("as_matrix", &vgl_rotation_3d<double>::as_matrix)
    .def("as_quaternion", &vgl_rotation_3d<double>::as_quaternion)
    .def("__repr__", streamToString<vgl_rotation_3d<double> >)
    .def(py::self * py::self);

  py::class_<vgl_pointset_3d<double> > (m, "pointset_3d")
    .def(py::init())
    .def(py::init<std::vector<vgl_point_3d<double> > >())
    .def(py::init<std::vector<vgl_point_3d<double> >, std::vector<vgl_vector_3d<double> > >())
    .def(py::init<std::vector<vgl_point_3d<double> >, std::vector<double> >())
    .def(py::init<std::vector<vgl_point_3d<double> >, std::vector<vgl_vector_3d<double> >,
        std::vector<double> >())
    .def("__len__", &vgl_pointset_3d<double>::size)
    .def("__repr__", [](vgl_pointset_3d<double> const& ptset){
        std::ostringstream buffer;
        buffer << std::boolalpha;
        buffer << "<vgl_pointset_3d";
        buffer << " n=" << ptset.size();
        buffer << " normals=" << ptset.has_normals();
        buffer << " scalars=" << ptset.has_scalars();
        buffer << ">";
        return buffer.str();
      })
    .def_property_readonly("has_normals", &vgl_pointset_3d<double>::has_normals)
    .def_property_readonly("has_scalars", &vgl_pointset_3d<double>::has_scalars)
    .def("add_point", &vgl_pointset_3d<double>::add_point)
    .def("add_point_with_normal", &vgl_pointset_3d<double>::add_point_with_normal)
    .def("add_point_with_scalar", &vgl_pointset_3d<double>::add_point_with_scalar)
    .def("add_point_with_normal_and_scalar", &vgl_pointset_3d<double>::add_point_with_normal_and_scalar)
    .def("points", &vgl_pointset_3d<double>::points)
    .def("normals", &vgl_pointset_3d<double>::normals)
    .def("scalars", &vgl_pointset_3d<double>::scalars)
    .def(py::self == py::self);

  py::class_<vgl_plane_3d<double> > (m, "plane_3d")
    .def(py::init())
    .def(py::init<double, double, double, double>())
    .def(py::init<vgl_vector_3d<double>, vgl_point_3d<double> >())
    .def(py::init<vgl_point_3d<double>, vgl_point_3d<double>, vgl_point_3d<double> >())
    .def("__repr__", streamToString<vgl_plane_3d<double> >)
    .def_property_readonly("a", &vgl_plane_3d<double>::a)
    .def_property_readonly("b", &vgl_plane_3d<double>::b)
    .def_property_readonly("c", &vgl_plane_3d<double>::c)
    .def_property_readonly("d", &vgl_plane_3d<double>::d)
    .def("set", &vgl_plane_3d<double>::set)
    .def_property_readonly("normal", &vgl_plane_3d<double>::normal)
    .def(py::self == py::self);

  py::class_<vgl_cylinder<double> > (m, "cylinder")
    .def(py::init())
    .def(py::init<vgl_point_3d<double>, double, double, vgl_vector_3d<double> >())
    .def("__repr__", streamToString<vgl_cylinder<double> >)
    .def_property("center", &vgl_cylinder<double>::center, &vgl_cylinder<double>::set_center)
    .def_property("radius", &vgl_cylinder<double>::radius, &vgl_cylinder<double>::set_radius)
    .def_property("length", &vgl_cylinder<double>::length, &vgl_cylinder<double>::set_length)
    .def_property("orientation", &vgl_cylinder<double>::orientation, &vgl_cylinder<double>::set_orientation)
    .def(py::self == py::self);

  py::class_<vgl_sphere_3d<double> > (m, "sphere_3d")
    .def(py::init())
    .def(py::init<vgl_point_3d<double>, double>())
    .def("__repr__", streamToString<vgl_sphere_3d<double> >)
    .def_property("center", &vgl_sphere_3d<double>::centre, &vgl_sphere_3d<double>::set_centre)
    .def_property("centre", &vgl_sphere_3d<double>::centre, &vgl_sphere_3d<double>::set_centre)
    .def_property("radius", &vgl_sphere_3d<double>::radius, &vgl_sphere_3d<double>::set_radius)
    .def(py::self == py::self);

  py::class_<vgl_polygon<double> > (m, "polygon")
    .def(py::init())
    .def(py::init<typename vgl_polygon<double>::sheet_t>())
    .def(py::init<std::vector<typename vgl_polygon<double>::sheet_t> >())
    .def("contains", (bool (vgl_polygon<double>::*)(vgl_point_2d<double> const&) const) &vgl_polygon<double>::contains)
    .def("__len__", &vgl_polygon<double>::num_sheets)
    .def("__getitem__", getitem_sheet<double>)
    .def("__repr__", [](vgl_polygon<double> const& p){
        std::ostringstream buffer;
        buffer << "<vgl_polygon num_sheets=" << p.num_sheets() << ">";
        return buffer.str();
    });

  py::class_<vgl_line_segment_2d<double> >(m, "line_segment_2d")
    .def(py::init())
    .def(py::init<vgl_point_2d<double>, vgl_point_2d<double> >())
    .def("__repr__", streamToString<vgl_line_segment_2d<double> >)
    .def_property_readonly("point1", &vgl_line_segment_2d<double>::point1)
    .def_property_readonly("point2", &vgl_line_segment_2d<double>::point2)
    .def("set", &vgl_line_segment_2d<double>::set)
    .def(py::self == py::self);

  py::class_<vgl_line_segment_3d<double> >(m, "line_segment_3d")
    .def(py::init())
    .def(py::init<vgl_point_3d<double>, vgl_point_3d<double> >())
    .def("__repr__", streamToString<vgl_line_segment_3d<double> >)
    .def_property_readonly("point1", &vgl_line_segment_3d<double>::point1)
    .def_property_readonly("point2", &vgl_line_segment_3d<double>::point2)
    .def("set", &vgl_line_segment_3d<double>::set)
    .def(py::self == py::self);

  py::class_<vgl_oriented_box_2d<double> >(m, "oriented_box_2d")
    .def(py::init())
    .def(py::init<vgl_line_segment_2d<double>, double>())
    .def("__repr__", streamToString<vgl_oriented_box_2d<double> >)
    .def_property_readonly("major_axis", &vgl_oriented_box_2d<double>::major_axis)
    .def_property_readonly("width", &vgl_oriented_box_2d<double>::width)
    .def_property_readonly("height", &vgl_oriented_box_2d<double>::height)
    .def("set", &vgl_oriented_box_2d<double>::set)
    .def(py::self == py::self);
}
}}
