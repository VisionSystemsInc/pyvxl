#include <pybind11/pybind11.h>

#ifdef PYVXL_WITH_CONTRIB_BPGL
  #include "pybpgl_algo.h"
#endif

namespace py = pybind11;

PYBIND11_MODULE(_vxl_contrib, m)
{
  m.doc() =  "Python bindings for the VXL computer vision libraries";

  py::module mod;

#ifdef PYVXL_WITH_CONTRIB_BPGL
  mod = m.def_submodule("bpgl");
  mod = mod.def_submodule("algo");
  pyvxl::bpgl::algo::wrap_bpgl_algo(mod);
#endif
}
