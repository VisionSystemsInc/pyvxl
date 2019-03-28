#include "pyvgl_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vgl/algo/vgl_compute_rigid_3d.h>
#include <vgl/algo/vgl_compute_similarity_3d.h>

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
