#include "pyvgl.h"
#include <vgl/vgl_vector_2d.h>
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_homg_point_2d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_homg_point_3d.h>
#include <vgl/vgl_ray_3d.h>
#include <vgl/vgl_pointset_3d.h>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_cylinder.h>
#include <vgl/vgl_sphere_3d.h>
#include <vgl/vgl_polygon.h>
#include <vgl/vgl_line_segment_2d.h>
#include <vgl/vgl_line_segment_3d.h>
#include <vgl/vgl_oriented_box_2d.h>
#include <vgl/vgl_fit_oriented_box_2d.h>
#include <vgl/vgl_box_2d.h>
#include <vgl/vgl_box_3d.h>
#include <vgl/vgl_intersection.h>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

// io classes for py::pickle
#include <vgl/io/vgl_io_point_2d.h>
#include <vgl/io/vgl_io_homg_point_2d.h>
#include <vgl/io/vgl_io_vector_2d.h>
#include <vgl/io/vgl_io_point_3d.h>
#include <vgl/io/vgl_io_homg_point_3d.h>
#include <vgl/io/vgl_io_vector_3d.h>
#include <vgl/io/vgl_io_polygon.h>

#include "../pyvxl_util.h"

#include <ios>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

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
double getitem_3d_homg(T const& a, long i)
{
  // wrap around
  if (i < 0) {
    i += 4;
  }
  if (i == 0) {
    return a.x();
  }
  else if (i == 1) {
    return a.y();
  }
  else if (i == 2) {
    return a.z();
  }
  else if (i == 3) {
    return a.w();
  }
  else {
    throw py::index_error("index out of range");
  }
  return 0;  // to avoid compiler warning
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

template<class T>
double getitem_2d_homg(T const& a, long i)
{
  // wrap around
  if (i < 0) {
    i += 3;
  }
  if (i == 0) {
    return a.x();
  }
  else if (i == 1) {
    return a.y();
  }
  else if (i == 2) {
    return a.w();
  }
  else {
    throw py::index_error("index out of range");
  }
  return 0;  // to avoid compiler warning
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

// Convert polygon to a list of sheets, where each sheet is a list of points, and each
// point is a list of T elements, instead of a vgl_point_2d vxl object -- this will
// auto-convert to python lists on the python side
template <typename T>
std::vector<std::vector<std::vector<T>>> vgl_polygon_as_array(vgl_polygon<T> const& p)
{
  // return an array of lists of [x,y] points
  std::vector<std::vector<std::vector<T>>> array_of_sheets;

  for (int i=0; i < p.num_sheets(); ++i) {

    // the vxl source sheet
    typename vgl_polygon<T>::sheet_t sheet = p[i];

    // create our target sheet
    std::vector<std::vector<T>> sheet_array;

    for (vgl_point_2d<T> p : sheet) {

      // create our target point
      std::vector<T> point_array;

      // Get the x and y values from the vgl_point_2d object
      point_array.push_back(p.x());
      point_array.push_back(p.y());

      // Add target point to target sheet
      sheet_array.push_back(point_array);
    }

    // Add target sheet to target return array
    array_of_sheets.push_back(sheet_array);
  }

  return array_of_sheets;
}

template<typename T>
void wrap_vgl_point_2d(py::module &m, std::string const& class_name)
{
  py::class_<vgl_point_2d<T> > (m, class_name.c_str())
    .def(py::init<>())
    .def(py::init<T,T>())
    .def(py::init(&type_from_buffer_2d<vgl_point_2d<T>, T>))
    .def("__len__", [](vgl_point_2d<T>){return (size_t)2;})
    .def("__getitem__", getitem_2d<vgl_point_2d<T> >)
    .def("__repr__", streamToString<vgl_point_2d<T> >)
    .def(py::pickle(&vslPickleGetState<vgl_point_2d<T> >,
                    &vslPickleSetState<vgl_point_2d<T> >))
    .def_property_readonly("x", (T (vgl_point_2d<T>::*)() const) &vgl_point_2d<T>::x)
    .def_property_readonly("y", (T (vgl_point_2d<T>::*)() const) &vgl_point_2d<T>::y)
    .def(py::self - py::self)
    .def(py::self == py::self)
    ;
}

template<typename T>
void wrap_vgl_homg_point_2d(py::module &m, std::string const& class_name)
{
  py::class_<vgl_homg_point_2d<T> > (m, class_name.c_str())
    .def(py::init<>())
    .def(py::init<T,T,T>(),
         py::arg("px"), py::arg("py"), py::arg("pw") = 1)
    .def("__len__", [](vgl_homg_point_2d<T>){return (size_t)3;})
    .def("__getitem__", getitem_2d_homg<vgl_homg_point_2d<T> >)
    .def("__repr__", streamToString<vgl_homg_point_2d<T> >)
    .def(py::pickle(&vslPickleGetState<vgl_homg_point_2d<T> >,
                    &vslPickleSetState<vgl_homg_point_2d<T> >))
    .def_property_readonly("x", (T (vgl_homg_point_2d<T>::*)() const) &vgl_homg_point_2d<T>::x)
    .def_property_readonly("y", (T (vgl_homg_point_2d<T>::*)() const) &vgl_homg_point_2d<T>::y)
    .def_property_readonly("w", (T (vgl_homg_point_2d<T>::*)() const) &vgl_homg_point_2d<T>::w)
    .def(py::self - py::self)
    .def(py::self == py::self)
    ;
}

template<typename T>
void wrap_vgl_vector_2d(py::module &m, std::string const& class_name)
{
  py::class_<vgl_vector_2d<T> > (m, class_name.c_str())
    .def(py::init<>())
    .def(py::init<T,T>())
    .def(py::init(&type_from_buffer_2d<vgl_vector_2d<T>, T>))
    .def("__len__", [](vgl_vector_2d<T>){return (size_t)2;})
    .def("__getitem__",getitem_2d<vgl_vector_2d<T> >)
    .def("__repr__", streamToString<vgl_vector_2d<T> >)
    .def(py::pickle(&vslPickleGetState<vgl_vector_2d<T> >,
                    &vslPickleSetState<vgl_vector_2d<T> >))
    .def_property_readonly("x", &vgl_vector_2d<T>::x)
    .def_property_readonly("y", &vgl_vector_2d<T>::y)
    .def("length", &vgl_vector_2d<T>::length)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self == py::self)
    .def("__neg__", [](vgl_vector_2d<T> const &v){return -v;})
    .def(py::self / T())
    .def(py::self * T())
    .def(T() * py::self)
    ;
}

template<typename T>
void wrap_vgl_point_3d(py::module &m, std::string const& class_name)
{
  py::class_<vgl_point_3d<T> > (m, class_name.c_str())
    .def(py::init<>())
    .def(py::init<T,T,T>())
    .def(py::init(&type_from_buffer_3d<vgl_point_3d<T>, T>))
    .def("__len__", [](vgl_point_3d<T>){return (size_t)3;})
    .def("__getitem__", getitem_3d<vgl_point_3d<T> >)
    .def("__repr__", streamToString<vgl_point_3d<T> >)
    .def(py::pickle(&vslPickleGetState<vgl_point_3d<T> >,
                    &vslPickleSetState<vgl_point_3d<T> >))
    .def_property_readonly("x", (T (vgl_point_3d<T>::*)() const) &vgl_point_3d<T>::x)
    .def_property_readonly("y", (T (vgl_point_3d<T>::*)() const) &vgl_point_3d<T>::y)
    .def_property_readonly("z", (T (vgl_point_3d<T>::*)() const) &vgl_point_3d<T>::z)
    .def(py::self - py::self)
    .def(py::self == py::self)
    ;
}

template<typename T>
void wrap_vgl_homg_point_3d(py::module &m, std::string const& class_name)
{
  py::class_<vgl_homg_point_3d<T> > (m, class_name.c_str())
    .def(py::init<>())
    .def(py::init<T,T,T,T>(),
         py::arg("px"), py::arg("py"), py::arg("pz"), py::arg("pw") = 1)
    .def("__len__", [](vgl_homg_point_3d<T>){return (size_t)4;})
    .def("__getitem__", getitem_3d_homg<vgl_homg_point_3d<T> >)
    .def("__repr__", streamToString<vgl_homg_point_3d<T> >)
    .def(py::pickle(&vslPickleGetState<vgl_homg_point_3d<T> >,
                    &vslPickleSetState<vgl_homg_point_3d<T> >))
    .def_property_readonly("x", (T (vgl_homg_point_3d<T>::*)() const) &vgl_homg_point_3d<T>::x)
    .def_property_readonly("y", (T (vgl_homg_point_3d<T>::*)() const) &vgl_homg_point_3d<T>::y)
    .def_property_readonly("z", (T (vgl_homg_point_3d<T>::*)() const) &vgl_homg_point_3d<T>::z)
    .def_property_readonly("w", (T (vgl_homg_point_3d<T>::*)() const) &vgl_homg_point_3d<T>::w)
    .def(py::self - py::self)
    .def(py::self == py::self)
    ;
}

template<typename T>
void wrap_vgl_vector_3d(py::module &m, std::string const& class_name)
{
  vgl_vector_3d<T> (*vector_3d_cross_product)(const vgl_vector_3d<T>&, const vgl_vector_3d<T>&) = &cross_product;
  py::class_<vgl_vector_3d<T> > (m, class_name.c_str())
    .def(py::init<>())
    .def(py::init<T,T,T>())
    .def(py::init(&type_from_buffer_3d<vgl_vector_3d<T>, T>))
    .def("__len__", [](vgl_vector_3d<T>){return (size_t)3;})
    .def("__getitem__", getitem_3d<vgl_vector_3d<T> >)
    .def("__repr__", streamToString<vgl_vector_3d<T> >)
    .def(py::pickle(&vslPickleGetState<vgl_vector_3d<T> >,
                    &vslPickleSetState<vgl_vector_3d<T> >))
    .def_property_readonly("x", &vgl_vector_3d<T>::x)
    .def_property_readonly("y", &vgl_vector_3d<T>::y)
    .def_property_readonly("z", &vgl_vector_3d<T>::z)
    .def("length", &vgl_vector_3d<T>::length)
    .def("cross_product", vector_3d_cross_product)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self == py::self)
    .def("__neg__", [](vgl_vector_3d<T> const &v){return -v;})
    .def(py::self / T())
    .def(py::self * T())
    .def(T() * py::self)
    ;
}

template<typename T>
void wrap_vgl_pointset_3d(py::module &m, std::string const& class_name)
{
    py::class_<vgl_pointset_3d<T> > (m, class_name.c_str())
    .def(py::init())
    .def(py::init<std::vector<vgl_point_3d<T> > >())
    .def(py::init<std::vector<vgl_point_3d<T> >, std::vector<vgl_vector_3d<T> > >())
    .def(py::init<std::vector<vgl_point_3d<T> >, std::vector<T> >())
    .def(py::init<std::vector<vgl_point_3d<T> >, std::vector<vgl_vector_3d<T> >,
        std::vector<T> >())
    .def("__len__", &vgl_pointset_3d<T>::size)
    .def("__repr__", [](vgl_pointset_3d<T> const& ptset){
        std::ostringstream buffer;
        buffer << std::boolalpha;
        buffer << "<vgl_pointset_3d";
        buffer << " n=" << ptset.size();
        buffer << " normals=" << ptset.has_normals();
        buffer << " scalars=" << ptset.has_scalars();
        buffer << ">";
        return buffer.str();
      })
    .def_property_readonly("has_normals", &vgl_pointset_3d<T>::has_normals)
    .def_property_readonly("has_scalars", &vgl_pointset_3d<T>::has_scalars)
    .def("add_point", &vgl_pointset_3d<T>::add_point)
    .def("add_point_with_normal", &vgl_pointset_3d<T>::add_point_with_normal)
    .def("add_point_with_scalar", &vgl_pointset_3d<T>::add_point_with_scalar)
    .def("add_point_with_normal_and_scalar", &vgl_pointset_3d<T>::add_point_with_normal_and_scalar)
    .def("points", &vgl_pointset_3d<T>::points)
    .def("normals", &vgl_pointset_3d<T>::normals)
    .def("scalars", &vgl_pointset_3d<T>::scalars)
    .def("append_pointset", &vgl_pointset_3d<T>::append_pointset)
    .def(py::self == py::self)
    .def("save", [](vgl_pointset_3d<T> const& ptset, std::string const& filename) {
         std::ofstream ofs(filename);
         if (!ofs.good()) {
           throw std::runtime_error("Bad filename");
         }
         ofs << ptset;
       })
    .def("load", [](vgl_pointset_3d<T> &ptset, std::string const& filename) {
         std::ifstream ifs(filename);
         if (!ifs.good()) {
           throw std::runtime_error("Bad filename");
         }
         ifs >> ptset;
       })
    ;
}

template<typename T>
void wrap_box_2d(py::module &m, std::string const& class_name)
{
  py::class_<vgl_box_2d<T> >(m, class_name.c_str())

    .def(py::init())
    // .def(py::init<double [], double []>())
    .def(py::init<vgl_point_2d<T>, vgl_point_2d<T> >())
    .def(py::init<T,T,T,T>(),
        py::arg("min_x"),py::arg("max_x"),py::arg("min_y"),py::arg("max_y"))
    .def("__repr__", streamToString<vgl_box_2d<T> >)

    .def_property("min_x", &vgl_box_2d<T>::min_x, &vgl_box_2d<T>::set_min_x)
    .def_property("min_y", &vgl_box_2d<T>::min_y, &vgl_box_2d<T>::set_min_y)
    .def_property("min_point", &vgl_box_2d<T>::min_point, &vgl_box_2d<T>::set_min_point)
    .def("set_min_position", &vgl_box_2d<T>::setmin_position)

    .def_property("max_x", &vgl_box_2d<T>::max_x, &vgl_box_2d<T>::set_max_x)
    .def_property("max_y", &vgl_box_2d<T>::max_y, &vgl_box_2d<T>::set_max_y)
    .def_property("max_point", &vgl_box_2d<T>::max_point, &vgl_box_2d<T>::set_max_point)
    .def("set_max_position", &vgl_box_2d<T>::setmax_position)

    .def_property("centroid_x", &vgl_box_2d<T>::centroid_x, &vgl_box_2d<T>::set_centroid_x)
    .def_property("centroid_y", &vgl_box_2d<T>::centroid_y, &vgl_box_2d<T>::set_centroid_y)
    .def_property("centroid", &vgl_box_2d<T>::centroid,
        (void (vgl_box_2d<T>::*)(vgl_point_2d<T> const&)) &vgl_box_2d<T>::set_centroid)
    .def("set_centroid_position", (void (vgl_box_2d<T>::*)(T const [])) &vgl_box_2d<T>::set_centroid)

    .def_property("width", &vgl_box_2d<T>::width, &vgl_box_2d<T>::set_width)
    .def_property("height", &vgl_box_2d<T>::height, &vgl_box_2d<T>::set_height)

    .def_property_readonly("volume", &vgl_box_2d<T>::volume)
    .def_property_readonly("is_empty", &vgl_box_2d<T>::is_empty)

    .def("add", (void (vgl_box_2d<T>::*)(vgl_point_2d<T> const&)) &vgl_box_2d<T>::add)
    .def("add", (void (vgl_box_2d<T>::*)(vgl_box_2d<T> const&)) &vgl_box_2d<T>::add)

    // .def_readonly("contains", (bool (vgl_box_2d<T>::*)(vgl_point_2d<T> const&) const) &vgl_box_2d<T>::contains)
    // .def_readonly("contains", (bool (vgl_box_2d<T>::*)(vgl_box_2d<T> const&) const) &vgl_box_2d<T>::contains)
    // .def_readonly("contains", (bool (vgl_box_2d<T>::*)(T const&, T const&) const) &vgl_box_2d<T>::contains)

    .def("expand_about_centroid", &vgl_box_2d<T>::expand_about_centroid)
    .def("scale_about_centroid", &vgl_box_2d<T>::scale_about_centroid)
    .def("scale_about_origin", &vgl_box_2d<T>::scale_about_origin)
    .def("empty", &vgl_box_2d<T>::empty)
    ;
}

template<typename T>
void wrap_box_3d(py::module &m, std::string const& class_name)
{
  py::class_<vgl_box_3d<T> >(m, class_name.c_str())

    .def(py::init())
    // .def(py::init<T [], T []>())
    .def(py::init<vgl_point_3d<T>, vgl_point_3d<T> >())
    .def(py::init<T,T,T,T,T,T>(),
        py::arg("min_x"),py::arg("min_y"),py::arg("min_z"),
        py::arg("max_x"),py::arg("max_y"),py::arg("max_z"))
    .def("__repr__", streamToString<vgl_box_3d<T> >)

    .def_property("min_x", &vgl_box_3d<T>::min_x, &vgl_box_3d<T>::set_min_x)
    .def_property("min_y", &vgl_box_3d<T>::min_y, &vgl_box_3d<T>::set_min_y)
    .def_property("min_z", &vgl_box_3d<T>::min_z, &vgl_box_3d<T>::set_min_z)
    .def_property("min_point", &vgl_box_3d<T>::min_point, &vgl_box_3d<T>::set_min_point)
    .def("set_min_position", &vgl_box_3d<T>::set_min_position)

    .def_property("max_x", &vgl_box_3d<T>::max_x, &vgl_box_3d<T>::set_max_x)
    .def_property("max_y", &vgl_box_3d<T>::max_y, &vgl_box_3d<T>::set_max_y)
    .def_property("max_z", &vgl_box_3d<T>::max_z, &vgl_box_3d<T>::set_max_z)
    .def_property("max_point", &vgl_box_3d<T>::max_point, &vgl_box_3d<T>::set_max_point)
    .def("set_max_position", &vgl_box_3d<T>::set_max_position)

    .def_property("centroid_x", &vgl_box_3d<T>::centroid_x, &vgl_box_3d<T>::set_centroid_x)
    .def_property("centroid_y", &vgl_box_3d<T>::centroid_y, &vgl_box_3d<T>::set_centroid_y)
    .def_property("centroid_z", &vgl_box_3d<T>::centroid_z, &vgl_box_3d<T>::set_centroid_z)
    .def_property("centroid", &vgl_box_3d<T>::centroid,
        (void (vgl_box_3d<T>::*)(vgl_point_3d<T> const&)) &vgl_box_3d<T>::set_centroid)
    .def("set_centroid_position", (void (vgl_box_3d<T>::*)(T const [])) &vgl_box_3d<T>::set_centroid)

    .def_property("width", &vgl_box_3d<T>::width, &vgl_box_3d<T>::set_width)
    .def_property("height", &vgl_box_3d<T>::height, &vgl_box_3d<T>::set_height)

    .def_property_readonly("volume", &vgl_box_3d<T>::volume)
    .def_property_readonly("vertices", &vgl_box_3d<T>::vertices)
    .def_property_readonly("is_empty", &vgl_box_3d<T>::is_empty)

    .def("add", (void (vgl_box_3d<T>::*)(vgl_point_3d<T> const&)) &vgl_box_3d<T>::add)
    .def("add", (void (vgl_box_3d<T>::*)(vgl_box_3d<T> const&)) &vgl_box_3d<T>::add)

    // .def_readonly("contains", (bool (vgl_box_3d<T>::*)(vgl_point_3d<T> const&) const) &vgl_box_3d<T>::contains)
    // .def_readonly("contains", (bool (vgl_box_3d<T>::*)(vgl_box_3d<T> const&) const) &vgl_box_3d<T>::contains)
    // .def_readonly("contains", (bool (vgl_box_3d<T>::*)(T const&, T const&, T const&) const) &vgl_box_3d<T>::contains)

    .def("expand_about_centroid", &vgl_box_3d<T>::expand_about_centroid)
    .def("scale_about_centroid", &vgl_box_3d<T>::scale_about_centroid)
    .def("scale_about_origin", &vgl_box_3d<T>::scale_about_origin)
    .def("empty", &vgl_box_3d<T>::empty)
    ;
}

void wrap_vgl(py::module &m)
{
  wrap_vgl_point_2d<double>(m, "point_2d");
  wrap_vgl_point_2d<float>(m, "point_2d_float");

  wrap_vgl_homg_point_2d<double>(m, "homg_point_2d");
  wrap_vgl_homg_point_2d<float>(m, "homg_point_2d_float");

  wrap_vgl_vector_2d<double>(m, "vector_2d");
  wrap_vgl_vector_2d<float>(m, "vector_2d_float");

  wrap_vgl_point_3d<double>(m, "point_3d");
  wrap_vgl_point_3d<float>(m, "point_3d_float");

  wrap_vgl_homg_point_3d<double>(m, "homg_point_3d");
  wrap_vgl_homg_point_3d<float>(m, "homg_point_3d_float");

  wrap_vgl_vector_3d<double>(m, "vector_3d");
  wrap_vgl_vector_3d<float>(m, "vector_3d_float");

  wrap_vgl_pointset_3d<double>(m, "pointset_3d");
  wrap_vgl_pointset_3d<float>(m, "pointset_3d_float");

  wrap_box_2d<double>(m,"box_2d");
  wrap_box_2d<float>(m,"box_2d_float");

  wrap_box_3d<double>(m,"box_3d");
  wrap_box_3d<float>(m,"box_3d_float");

  py::class_<vgl_ray_3d<double> >(m, "ray_3d")
    .def(py::init<vgl_point_3d<double>, vgl_vector_3d<double> >())
    .def_property_readonly("origin", &vgl_ray_3d<double>::origin)
    .def_property_readonly("direction", &vgl_ray_3d<double>::direction);


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
    .def("as_array", &vgl_polygon_as_array<double>)
    .def("__len__", &vgl_polygon<double>::num_sheets)
    .def("__getitem__", getitem_sheet<double>)
    .def("push_back_xy",
        overload_cast_<double, double>()(&vgl_polygon<double>::push_back))
    .def("push_back_point",
        overload_cast_<vgl_point_2d<double> const&>()(&vgl_polygon<double>::push_back))
    .def("push_back_sheet",
        overload_cast_<std::vector<vgl_point_2d<double> > const&>()(&vgl_polygon<double>::push_back))
    .def(py::pickle(&vslPickleGetState<vgl_polygon<double>>,
                    &vslPickleSetState<vgl_polygon<double>>))
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

  vgl_oriented_box_2d<double> (vgl_fit_oriented_box_2d<double>::*oriented_fit_box)() = &vgl_fit_oriented_box_2d<double>::fitted_box;
  py::class_<vgl_fit_oriented_box_2d<double> >(m, "fit_oriented_box_2d")
    .def(py::init<vgl_polygon<double> >())
    .def("fitted_box", oriented_fit_box);

  m.def("intersection", [](vgl_ray_3d<double> const& r, vgl_plane_3d<double> const& p) {
        vgl_point_3d<double> pt;
        if(!vgl_intersection(r,p,pt)) {
          throw std::runtime_error("ray does not intersect plane");
        }
        return pt;
        });
}
}}

PYBIND11_MODULE(_vgl, m)
{
  m.doc() =  "Python bindings for the VGL computer vision libraries";

  pyvxl::vgl::wrap_vgl(m);


  /* py::module mod = m.def_submodule("vgl"); */
  /* pyvxl::vgl::wrap_vgl(mod); */
  /* mod = mod.def_submodule("algo"); */
  /* pyvxl::vgl::algo::wrap_vgl_algo(mod); */
}
