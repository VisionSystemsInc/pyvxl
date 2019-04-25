#include "pysdet_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <string>
#include <iostream>

#include <sdet/sdet_texture_classifier.h>
#include <sdet/sdet_texture_classifier_params.h>
#include <sdet/algo/sdet_classify_clouds.h>

namespace py = pybind11;

namespace pyvxl { namespace sdet { namespace algo {

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

void wrap_sdet_algo(py::module &m)
{
  m.def("load_classifier", &load_classifier,
        py::arg("input_ins_path"));

  m.def("classify_clouds", &sdet_classify_clouds,
        py::arg("cloud_classifier"), py::arg("texton_dict_path"),
        py::arg("image"), py::arg("i"), py::arg("j"),
        py::arg("ni"), py::arg("nj"), py::arg("block_size"),
        py::arg("cat_ids_file"), py::arg("first_category"),
        py::arg("scale_factor"));
}

}}}

PYBIND11_MODULE(_sdet_algo, m)
{
  m.doc() = "Python bindings for the VXL SDET algo computer vision library";

  pyvxl::sdet::algo::wrap_sdet_algo(m);
}

