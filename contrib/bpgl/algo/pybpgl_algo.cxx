#include "pybpgl_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>
#include <vector>
#include <array>

#include <vpgl/vpgl_affine_camera.h>
#include <bpgl/algo/bpgl_heightmap_from_disparity.h>
#include <bpgl/algo/bpgl_rectify_affine_image_pair.h>

namespace py = pybind11;

namespace pyvxl {
namespace bpgl {
namespace algo {
  void wrap_bpgl_algo(py::module &m)
  {
    m.def("heightmap_from_disparity_affine", &bpgl_heightmap_from_disparity<vpgl_affine_camera<double> >);
    py::class_<bpgl_rectify_affine_image_pair>(m, "bpgl_rectify_affine_image_pair")
      .def(py::init<>())
      .def(py::init<vil_image_view_base_sptr const&, vpgl_affine_camera<double> const&, vil_image_view_base_sptr const&, vpgl_affine_camera<double> const&>())
      .def(py::init<vil_image_view<unsigned char> const&, vpgl_affine_camera<double> const&, vil_image_view<unsigned char> const&, vpgl_affine_camera<double> const&>())
      .def("set_param_values", &bpgl_rectify_affine_image_pair::set_param_values, "set all parameters to avoid wrapping a separate param struct",
          py::arg("min_disparity_z") = NAN, py::arg("n_points") = 1000, py::arg("upsample_scale") = 1.0)
      .def("process", &bpgl_rectify_affine_image_pair::process, "main process method")
      .def("rectified_fview0", &bpgl_rectify_affine_image_pair::rectified_fview0, "rectified image 0 of a pair")
      .def("rectified_fview1", &bpgl_rectify_affine_image_pair::rectified_fview1, "rectified image 1 of a pair");
  }}}
}

PYBIND11_MODULE(_bpgl_algo, m)
{
  m.doc() =  "Python bindings for the VXL computer vision libraries";

  pyvxl::bpgl::algo::wrap_bpgl_algo(m);
}

