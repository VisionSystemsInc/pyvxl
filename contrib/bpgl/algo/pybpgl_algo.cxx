#include "pybpgl_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>
#include <vector>
#include <array>

#include <vil/vil_image_view.h>
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
    
    py::class_<bpgl_rectify_affine_image_pair> (m, "bpgl_rectify_affine_image_pair")
      .def(py::init<>())
      .def(py::init<vil_image_view<unsigned char> const&, 
                    vpgl_affine_camera<double> const&, 
                    vil_image_view<unsigned char> const&, 
                    vpgl_affine_camera<double> const&>())
      .def("set_param_values", &bpgl_rectify_affine_image_pair::set_param_values)
      .def("process", &bpgl_rectify_affine_image_pair::process)
      .def_property_readonly("acam0", &bpgl_rectify_affine_image_pair::acam0)
      .def_property_readonly("acam1", &bpgl_rectify_affine_image_pair::acam1)
      .def_property_readonly("rect_acam0", &bpgl_rectify_affine_image_pair::rect_acam0)
      .def_property_readonly("rect_acam1", &bpgl_rectify_affine_image_pair::rect_acam1)
      .def_property_readonly("input_float_view0", &bpgl_rectify_affine_image_pair::input_float_view0)
      .def_property_readonly("input_float_view1", &bpgl_rectify_affine_image_pair::input_float_view1)
      .def_property_readonly("rectified_fview0", &bpgl_rectify_affine_image_pair::rectified_fview0)
      .def_property_readonly("rectified_fview1", &bpgl_rectify_affine_image_pair::rectified_fview1)
      .def_property_readonly("H0", &bpgl_rectify_affine_image_pair::H0)
      .def_property_readonly("H1", &bpgl_rectify_affine_image_pair::H1)

      ;

  }
}}}

PYBIND11_MODULE(_bpgl_algo, m)
{
  m.doc() =  "Python bindings for the VXL computer vision libraries";

  pyvxl::bpgl::algo::wrap_bpgl_algo(m);
}

