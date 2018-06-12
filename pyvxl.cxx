#include <pybind11/pybind11.h>
#include "pyvnl.h"
#include "pyvgl.h"
#include "pyvpgl.h"
#include "pyvil.h"
#include "pyvgl_algo.h"
#include "pyvpgl_algo.h"

namespace py = pybind11;

PYBIND11_MODULE(vxl, m)
{
  m.doc() =  "Python bindings for the VXL computer vision libraries";

  py::module mod = m.def_submodule("vnl");
  pyvxl::vnl::wrap_vnl(mod);

  mod = m.def_submodule("vgl");
  pyvxl::vgl::wrap_vgl(mod);
  mod = mod.def_submodule("algo");
  pyvxl::vgl::algo::wrap_vgl_algo(mod);

  mod = m.def_submodule("vpgl");
  pyvxl::wrap_vpgl(mod);
  mod = mod.def_submodule("algo");
  pyvxl::vpgl::algo::wrap_vpgl_algo(mod);

  mod = m.def_submodule("vil");
  pyvxl::vil::wrap_vil(mod);
}
