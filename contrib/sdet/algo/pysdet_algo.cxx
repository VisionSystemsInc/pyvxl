#include "pysdet_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <sdet/algo/sdet_classify.h>

#include "pyvxl_holder_types.h"

namespace py = pybind11;

namespace pyvxl { namespace sdet { namespace algo {

void wrap_sdet_algo(py::module &m)
{

  m.def("classify", &sdet_classify, py::call_guard<py::gil_scoped_release>(),
        py::arg("classifier"),
        py::arg("float_image"),
        py::arg("category"),
        py::arg("cat_ids_file") = std::string(""),
        py::arg("scale_factor") = 1.0 / 2048.0);
}

}}}

PYBIND11_MODULE(_sdet_algo, m)
{
  m.doc() = "Python bindings for the VXL SDET algo computer vision library";

  pyvxl::sdet::algo::wrap_sdet_algo(m);
}

