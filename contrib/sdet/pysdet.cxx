#include "pysdet.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <string>
#include <iostream>

#include <sdet/sdet_atmospheric_image_classifier.h>
#include <sdet/sdet_texture_classifier.h>
#include <sdet/sdet_texture_classifier_params.h>

namespace py = pybind11;

namespace pyvxl { namespace sdet {

sdet_texture_classifier load_classifier(std::string input_ins_path)
{

  sdet_texture_classifier_params dummy;
  sdet_texture_classifier tc(dummy);

  tc.load_data(input_ins_path);
  std::cout << " loaded classifier with params: " << tc << std::endl;

  tc.filter_responses().set_params(tc.n_scales_,tc.scale_interval_,tc.lambda0_,tc.lambda1_,tc.angle_interval_,tc.cutoff_per_);

  std::cout << " in the loaded classifier max filter radius: " << tc.max_filter_radius() << std::endl;

  return tc;
}

sdet_atmospheric_image_classifier load_atmospheric_classifier(std::string input_ins_path)
{
  sdet_texture_classifier tc = load_classifier(input_ins_path);
  return sdet_atmospheric_image_classifier(tc);
}

void wrap_sdet(py::module &m)
{
  py::class_<sdet_texture_classifier> (m, "texture_classifier")
    .def(py::init<sdet_texture_classifier_params const&>())
    .def("load_dictionary", &sdet_texture_classifier::load_dictionary)
    .def("max_filter_radius", &sdet_texture_classifier::max_filter_radius)
    .def_readwrite("block_size", &sdet_texture_classifier::block_size_);

  py::class_<sdet_atmospheric_image_classifier, sdet_texture_classifier /* <- Parent */ > (m, "atmospheric_image_classifier")
    .def(py::init<sdet_texture_classifier_params const&>());


  m.def("load_classifier", &load_classifier, py::arg("input_ins_path"));
  m.def("load_atmospheric_classifier", &load_atmospheric_classifier, py::arg("input_ins_path"));
}

}}

PYBIND11_MODULE(_sdet, m)
{
  m.doc() = "Python bindings for the VXL SDET computer vision library";

  pyvxl::sdet::wrap_sdet(m);
}

