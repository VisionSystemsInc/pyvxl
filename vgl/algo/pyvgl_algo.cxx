#include "pyvgl_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vgl/algo/vgl_compute_rigid_3d.h>
#include <vgl/algo/vgl_compute_similarity_3d.h>
#include <vgl/algo/vgl_rotation_3d.h>

#include "../../pyvxl_util.h"

namespace py = pybind11;

namespace pyvxl { namespace vgl { namespace algo {

void wrap_vgl_algo(py::module &m)
{

  py::class_ <vgl_compute_similarity_3d<double> > (m, "compute_similarity_3d")
    .def(py::init<std::vector<vgl_point_3d<double> >, std::vector<vgl_point_3d<double> > >())
    .def("estimate", &vgl_compute_similarity_3d<double>::estimate)
    .def("rotation", &vgl_compute_similarity_3d<double>::rotation)
    .def("translation", &vgl_compute_similarity_3d<double>::translation)
    .def("scale", &vgl_compute_similarity_3d<double>::scale);

  py::class_ <vgl_compute_rigid_3d<double> > (m, "compute_rigid_3d")
    .def(py::init<std::vector<vgl_point_3d<double> >, std::vector<vgl_point_3d<double> > >())
    .def("estimate", &vgl_compute_rigid_3d<double>::estimate)
    .def("rotation", &vgl_compute_rigid_3d<double>::rotation)
    .def("translation", &vgl_compute_rigid_3d<double>::translation);

  py::class_ <vgl_rotation_3d<double> > (m, "rotation_3d")
    .def(py::init<>())
    .def(py::init<vnl_vector_fixed<double,4> >())
    .def(py::init<vnl_matrix_fixed<double,3,3> >())
    .def(py::init<vnl_vector_fixed<double,3> >())
    .def(py::init<vnl_vector_fixed<double,3>,vnl_vector_fixed<double,3> >())
    .def("as_matrix", &vgl_rotation_3d<double>::as_matrix)
    .def("as_quaternion", &vgl_rotation_3d<double>::as_quaternion)
    .def("inverse", &vgl_rotation_3d<double>::inverse)
    .def("transpose", &vgl_rotation_3d<double>::transpose)
    .def("__repr__", streamToString<vgl_rotation_3d<double> >)
    .def(py::self * vgl_vector_3d<double>())
    .def(py::self * vgl_point_3d<double>())
    .def(py::self * py::self)
    .def(py::self == py::self)
    ;

}
}}}

PYBIND11_MODULE(_vgl_algo, m)
{
  m.doc() =  "Python bindings for the VGL Algo computer vision libraries";

  pyvxl::vgl::algo::wrap_vgl_algo(m);


  /* py::module mod = m.def_submodule("vgl"); */
  /* pyvxl::vgl::wrap_vgl(mod); */
  /* mod = mod.def_submodule("algo"); */
  /* pyvxl::vgl::algo::wrap_vgl_algo(mod); */
}
