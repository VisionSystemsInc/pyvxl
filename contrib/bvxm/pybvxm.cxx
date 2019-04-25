#include "pybvxm.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <bvxm/bvxm_edge_util.h>

namespace py = pybind11;

namespace pyvxl { namespace bvxm {

void wrap_bvxm(py::module &m)
{
  m.def("detect_edges", &bvxm_edge_util::detect_edges,
        py::arg("input_img"), py::arg("noise_multiplier"), py::arg("smooth"),
        py::arg("automatic_threshold"), py::arg("junctionp"), py::arg("aggressive_junction_closure"));
}

}}

PYBIND11_MODULE(_bvxm, m)
{
  m.doc() = "Python bindings for the VXL BVXM computer vision library";

  pyvxl::bvxm::wrap_bvxm(m);
}

