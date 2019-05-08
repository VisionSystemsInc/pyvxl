#include "pysdet_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <sdet/algo/sdet_classify_clouds.h>

#include "pyvxl_holder_types.h"

namespace py = pybind11;

namespace pyvxl { namespace sdet { namespace algo {

void wrap_sdet_algo(py::module &m)
{
  m.def("classify_clouds", &sdet_classify_clouds,
        py::arg("cloud_classifier"), py::arg("texton_dict_path"),
        py::arg("image_resource"),
        py::arg("i"), py::arg("j"),
        py::arg("ni"), py::arg("nj"),
        py::arg("block_size"),
        py::arg("first_category"),
        py::arg("cat_ids_file") = std::string(""),
        py::arg("scale_factor") = 1.0 / 2048.0);
}

}}}

PYBIND11_MODULE(_sdet_algo, m)
{
  m.doc() = "Python bindings for the VXL SDET algo computer vision library";

  pyvxl::sdet::algo::wrap_sdet_algo(m);
}

