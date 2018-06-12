#include "pyvgl.h"
#include <vgl/vgl_vector_2d.h>
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_point_3d.h>
#include <vgl/algo/vgl_rotation_3d.h>

#include <vgl/vgl_homg_point_2d.h>
#include <vgl/vgl_homg_point_3d.h>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_homg_line_3d_2_points.h>
#include <vgl/algo/vgl_homg_operators_3d.h>
#include <vgl/vgl_homg_plane_3d.h>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vgl/algo/vgl_h_matrix_2d.h>
#include <vgl/vgl_closest_point.h>
#include <vgl/vgl_distance.h>
#include <vgl/vgl_ray_3d.h>

namespace py = pybind11;

namespace pyvxl {
    namespace vgl {

        // Helper functions

        template<class T>
        double getitem_4d(T const &a, long i) {
          // wrap around
          if (i < 0) {
            i += 3;
          }
          if (i == 0) {
            return a.a();
          } else if (i == 1) {
            return a.b();
          } else if (i == 2) {
            return a.c();
          } else if (i == 3) {
            return a.d();
          } else {
            throw py::index_error("index out of range");
          }
          return 0; // to avoid compiler warning
        }

        template<class T>
        double getitem_3d(T const &a, long i) {
          // wrap around
          if (i < 0) {
            i += 3;
          }
          if (i == 0) {
            return a.x();
          } else if (i == 1) {
            return a.y();
          } else if (i == 2) {
            return a.z();
          } else {
            throw py::index_error("index out of range");
          }
          return 0; // to avoid compiler warning
        }

        template<class T>
        double getitem_2d(T const &a, long i) {
          // wrap around
          if (i < 0) {
            i += 2;
          }
          if (i == 0) {
            return a.x();
          } else if (i == 1) {
            return a.y();
          } else {
            throw py::index_error("index out of range");
          }
          return 0; // to avoid compiler warning
        }

        template<class T>
        double getitem_homg_2d(T const &a, long i) {
          // wrap around
          if (i < 0) {
            i += 2;
          }
          if (i == 0) {
            return a.x();
          } else if (i == 1) {
            return a.y();
          } else if (i == 2) {
            return a.w();
          } else {
            throw py::index_error("index out of range");
          }
          return 0; // to avoid compiler warning
        }

        template<class T>
        double getitem_homg_3d(T const &a, long i) {
          // wrap around
          if (i < 0) {
            i += 3;
          }
          if (i == 0) {
            return a.x();
          } else if (i == 1) {
            return a.y();
          } else if (i == 2) {
            return a.z();
          } else if (i == 3) {
            return a.w();
          } else {
            throw py::index_error("index out of range");
          }
          return 0; // to avoid compiler warning
        }

        template<class T, class BUFF_T>
        T type_from_buffer_2d(py::array_t<BUFF_T> b) {
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
          const BUFF_T *data_ptr = static_cast<BUFF_T *>(info.ptr);
          const size_t stride = info.strides[0] / sizeof(BUFF_T);
          BUFF_T x = *data_ptr;
          BUFF_T y = *(data_ptr + stride);

          return T(x, y);
        }

        template<class T, class BUFF_T>
        T type_from_buffer_3d(py::array_t<BUFF_T> b) {
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
          const BUFF_T *data_ptr = static_cast<BUFF_T *>(info.ptr);
          const size_t stride = info.strides[0] / sizeof(BUFF_T);
          BUFF_T x = *data_ptr;
          BUFF_T y = *(data_ptr + stride);
          BUFF_T z = *(data_ptr + 2 * stride);

          return T(x, y, z);
        }

        template<class T, class BUFF_T>
        T type_from_buffer_4d(py::array_t<BUFF_T> buffer) {
          py::buffer_info info = buffer.request();
          if (info.format != py::format_descriptor<BUFF_T>::format()) {
            throw std::runtime_error("Incompatible scalar type");
          }
          if (info.ndim != 1) {
            throw std::runtime_error("Expecting a 1-dimensional vector");
          }
          const size_t num_elements = info.shape[0];
          if (num_elements != 4) {
            throw std::runtime_error("Expecting 4-d input vector");
          }
          // in-place constructor
          const BUFF_T *data_ptr = static_cast<BUFF_T *>(info.ptr);
          const size_t stride = info.strides[0] / sizeof(BUFF_T);
          BUFF_T a = *data_ptr;
          BUFF_T b = *(data_ptr + stride);
          BUFF_T c = *(data_ptr + 2 * stride);
          BUFF_T d = *(data_ptr + 3 * stride);

          return T(a, b, c, d);
        }

        template<class T, class BUFF_T>
        T type_from_buffer_homg_2d(py::array_t<BUFF_T> b) {
          py::buffer_info info = b.request();
          if (info.format != py::format_descriptor<BUFF_T>::format()) {
            throw std::runtime_error("Incompatible scalar type");
          }
          if (info.ndim != 1) {
            throw std::runtime_error("Expecting a 1-dimensional vector");
          }
          const size_t num_elements = info.shape[0];
          if (num_elements > 3) {
            throw std::runtime_error("Expecting 2-d or 3-d input vector");
          }
          // in-place constructor
          const BUFF_T *data_ptr = static_cast<BUFF_T *>(info.ptr);
          const size_t stride = info.strides[0] / sizeof(BUFF_T);
          BUFF_T x = *data_ptr;
          BUFF_T y = *(data_ptr + stride);

          if (num_elements == 2)
            return T(x, y);

          BUFF_T w = *(data_ptr + 2 * stride);

          return T(x, y, w);
        }

        template<class T, class BUFF_T>
        T type_from_buffer_homg_3d(py::array_t<BUFF_T> b) {
          py::buffer_info info = b.request();
          if (info.format != py::format_descriptor<BUFF_T>::format()) {
            throw std::runtime_error("Incompatible scalar type");
          }
          if (info.ndim != 1) {
            throw std::runtime_error("Expecting a 1-dimensional vector");
          }
          const size_t num_elements = info.shape[0];
          if (num_elements == 2 || num_elements > 4) {
            throw std::runtime_error("Expecting 3-d or 4-d input vector");
          }
          // in-place constructor
          const BUFF_T *data_ptr = static_cast<BUFF_T *>(info.ptr);
          const size_t stride = info.strides[0] / sizeof(BUFF_T);
          BUFF_T x = *data_ptr;
          BUFF_T y = *(data_ptr + stride);
          BUFF_T z = *(data_ptr + 2 * stride);

          if (num_elements == 3)
            return T(x, y, z);

          BUFF_T w = *(data_ptr + 3 * stride);

          return T(x, y, z, w);
        }

        void wrap_vgl(py::module &m) {
          py::class_<vgl_point_2d<double> >(m, "point_2d")
              .def(py::init<double, double>())
              .def(py::init(&type_from_buffer_2d<vgl_point_2d<double>, double>))
              .def("__len__", [](vgl_point_2d<double>) { return (size_t) 2; })
              .def("__getitem__", getitem_2d<vgl_point_2d<double> >)
              .def_property_readonly("x", (double (vgl_point_2d<double>::*)() const) &vgl_point_2d<double>::x)
              .def_property_readonly("y", (double (vgl_point_2d<double>::*)() const) &vgl_point_2d<double>::y)
              .def(py::self - py::self);

          py::class_<vgl_vector_2d<double> >(m, "vector_2d")
              .def(py::init<double, double>())
              .def(py::init(&type_from_buffer_2d<vgl_vector_2d<double>, double>))
              .def("__len__", [](vgl_vector_2d<double>) { return (size_t) 2; })
              .def("__getitem__", getitem_2d<vgl_vector_2d<double> >)
              .def_property_readonly("x", &vgl_vector_2d<double>::x)
              .def_property_readonly("y", &vgl_vector_2d<double>::y)
              .def("length", &vgl_vector_2d<double>::length)
              .def(py::self + py::self)
              .def(py::self - py::self);

          py::class_<vgl_h_matrix_2d<double> >(m, "matrix_2d")
              .def(py::init<>())
              .def(py::self * vgl_homg_point_2d<double>())
              .def("get_matrix", &vgl_h_matrix_2d<double>::get_matrix);

          py::class_<vgl_point_3d<double> >(m, "point_3d")
              .def(py::init<double, double, double>())
              .def(py::init(&type_from_buffer_3d<vgl_point_3d<double>, double>))
              .def("__len__", [](vgl_point_3d<double>) { return (size_t) 3; })
              .def("__getitem__", getitem_3d<vgl_point_3d<double> >)
              .def_property_readonly("x", (double (vgl_point_3d<double>::*)() const) &vgl_point_3d<double>::x)
              .def_property_readonly("y", (double (vgl_point_3d<double>::*)() const) &vgl_point_3d<double>::y)
              .def_property_readonly("z", (double (vgl_point_3d<double>::*)() const) &vgl_point_3d<double>::z)
              .def(py::self - py::self);

          py::class_<vgl_vector_3d<double> >(m, "vector_3d")
              .def(py::init<double, double, double>())
              .def(py::init(&type_from_buffer_3d<vgl_vector_3d<double>, double>))
              .def("__len__", [](vgl_vector_3d<double>) { return (size_t) 3; })
              .def("__getitem__", getitem_3d<vgl_vector_3d<double> >)
              .def_property_readonly("x", &vgl_vector_3d<double>::x)
              .def_property_readonly("y", &vgl_vector_3d<double>::y)
              .def_property_readonly("z", &vgl_vector_3d<double>::z)
              .def("length", &vgl_vector_3d<double>::length)
              .def(py::self + py::self)
              .def(py::self - py::self);

          py::class_<vgl_rotation_3d<double> >(m, "rotation_3d")
              .def(py::init<vnl_vector_fixed<double, 4> >())
              .def(py::init<vnl_matrix_fixed<double, 3, 3> >())
              .def("as_matrix", &vgl_rotation_3d<double>::as_matrix)
              .def("as_quaternion", &vgl_rotation_3d<double>::as_quaternion)
              .def(py::self * py::self);

          py::class_<vgl_homg_point_2d<double> >(m, "homg_point_2d")
              .def(py::init<double, double>())
              .def(py::init<double, double, double>())
              .def(py::init(&type_from_buffer_homg_2d<vgl_homg_point_2d<double>, double>))
              .def("__len__", [](vgl_homg_point_2d<double>) { return (size_t) 3; })
              .def("__getitem__", getitem_homg_2d<vgl_homg_point_2d<double> >)
              .def_property_readonly("x", (double (vgl_homg_point_2d<double>::*)() const) &vgl_homg_point_2d<double>::x)
              .def_property_readonly("y", (double (vgl_homg_point_2d<double>::*)() const) &vgl_homg_point_2d<double>::y)
              .def_property_readonly("w", (double (vgl_homg_point_2d<double>::*)() const) &vgl_homg_point_2d<double>::w)
              .def(py::self - py::self);

          py::class_<vgl_homg_point_3d<double> >(m, "homg_point_3d")
              .def(py::init<double, double, double>())
              .def(py::init<double, double, double, double>())
              .def(py::init(&type_from_buffer_homg_3d<vgl_homg_point_3d<double>, double>))
              .def("__len__", [](vgl_homg_point_3d<double>) { return (size_t) 4; })
              .def("__getitem__", getitem_homg_3d<vgl_homg_point_3d<double> >)
              .def_property_readonly("x", (double (vgl_homg_point_3d<double>::*)() const) &vgl_homg_point_3d<double>::x)
              .def_property_readonly("y", (double (vgl_homg_point_3d<double>::*)() const) &vgl_homg_point_3d<double>::y)
              .def_property_readonly("z", (double (vgl_homg_point_3d<double>::*)() const) &vgl_homg_point_3d<double>::z)
              .def_property_readonly("w", (double (vgl_homg_point_3d<double>::*)() const) &vgl_homg_point_3d<double>::w)
              .def(py::self - py::self);

          py::class_<vgl_plane_3d<double> >(m, "plane_3d")
              .def(py::init<double, double, double, double>())
              .def(py::init(&type_from_buffer_4d<vgl_plane_3d<double>, double>))
              .def(py::init<const vgl_vector_3d<double> &, const vgl_point_3d<double> &>())
              .def("__len__", [](vgl_plane_3d<double>) { return (size_t) 4; })
              .def("__getitem__", getitem_4d<vgl_plane_3d<double> >)
              .def_property_readonly("a", (double (vgl_plane_3d<double>::*)() const) &vgl_plane_3d<double>::a)
              .def_property_readonly("b", (double (vgl_plane_3d<double>::*)() const) &vgl_plane_3d<double>::b)
              .def_property_readonly("c", (double (vgl_plane_3d<double>::*)() const) &vgl_plane_3d<double>::c)
              .def_property_readonly("d", (double (vgl_plane_3d<double>::*)() const) &vgl_plane_3d<double>::d);

          py::class_ <vgl_homg_plane_3d<double> >(m, "homg_plane_3d")
              .def(py::init<>())
              .def(py::init<double, double, double, double>())
              .def(py::init(&type_from_buffer_4d<vgl_homg_plane_3d<double>, double>))
              .def(py::init<const vgl_homg_plane_3d<double> &>())
              .def("__len__", [](vgl_homg_plane_3d<double>) { return (size_t) 4; })
              .def("__getitem__", getitem_4d<vgl_plane_3d<double> >)
              .def_property_readonly("a", (double (vgl_homg_plane_3d<double>::*)() const) &vgl_plane_3d<double>::a)
              .def_property_readonly("b", (double (vgl_homg_plane_3d<double>::*)() const) &vgl_plane_3d<double>::b)
              .def_property_readonly("c", (double (vgl_homg_plane_3d<double>::*)() const) &vgl_plane_3d<double>::c)
              .def_property_readonly("d", (double (vgl_homg_plane_3d<double>::*)() const) &vgl_plane_3d<double>::d)
              .def_property_readonly("nx", (double (vgl_homg_plane_3d<double>::*)() const) &vgl_plane_3d<double>::nx)
              .def_property_readonly("ny", (double (vgl_homg_plane_3d<double>::*)() const) &vgl_plane_3d<double>::ny)
              .def_property_readonly("nz", (double (vgl_homg_plane_3d<double>::*)() const) &vgl_plane_3d<double>::nz);

          py::class_<vgl_homg_line_3d_2_points<double> >(m, "homg_line_3d_2_points")
              .def(py::init<>())
              .def(py::init<vgl_homg_point_3d<double>, vgl_homg_point_3d<double> >())
              .def_property_readonly("point_finite", &vgl_homg_line_3d_2_points<double>::point_finite)
              .def_property_readonly("point_infinite", &vgl_homg_line_3d_2_points<double>::point_infinite);

          py::class_<vgl_homg_operators_3d<double> >(m, "homg_operators_3d")
              .def(py::init<>())
              .def_static("intersect_line_and_plane", &vgl_homg_operators_3d<double>::intersect_line_and_plane)
              .def_static("perp_projection", &vgl_homg_operators_3d<double>::perp_projection);

          py::class_<vgl_line_3d_2_points<double> >(m, "line_3d_2_points")
              .def(py::init<>())
              .def(py::init<vgl_point_3d<double>, vgl_point_3d<double> >())
              .def_property_readonly("point1", &vgl_line_3d_2_points<double>::point1)
              .def_property_readonly("point2", &vgl_line_3d_2_points<double>::point2);

          py::class_<vgl_ray_3d<double> >(m, "ray_3d")
              .def(py::init<>())
              .def_property_readonly("origin", &vgl_ray_3d<double>::origin)
              .def_property_readonly("direction", &vgl_ray_3d<double>::direction);

          m.def("closest_point", (vgl_point_3d<double>
              (*)(vgl_line_3d_2_points<double> const&, vgl_point_3d<double> const&)) &vgl_closest_point<double>);

          m.def("distance", (double
          (*)(vgl_line_3d_2_points<double> const&, vgl_point_3d<double> const&)) &vgl_distance<double>);

          m.def("distance", (double
          (*)(vgl_point_3d<double> const&, vgl_point_3d<double> const&)) &vgl_distance<double>);
        }
    }
}
