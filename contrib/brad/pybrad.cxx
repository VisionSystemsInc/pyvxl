#include "pybrad.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>
#include <vector>
#include <array>

#include <vil/vil_image_view.h>
#include <brad/brad_image_atmospherics_est.h>

namespace py = pybind11;

namespace pyvxl { namespace brad {

void wrap_brad(py::module &m)
{
  m.def("estimate_reflectance", &brad_estimate_reflectance_image);
}

}}

PYBIND11_MODULE(_brad, m)
{
  m.doc() =  "Python bindings for the VXL BRAD computer vision library";

  pyvxl::brad::wrap_brad(m);

}
