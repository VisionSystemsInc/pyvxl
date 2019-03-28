#include <pybind11/pybind11.h>
#include "pyvnl.h"
#include "pyvgl.h"
#include "pyvpgl.h"
#include "pyvil.h"
#include "pyvgl_algo.h"
#include "pyvpgl_algo.h"

#ifdef PYVXL_WITH_CONTRIB_BPGL
#include "pybpgl_algo.h"
#endif

namespace py = pybind11;

// helper function to check if py::module import exists
bool import_exists(std::string const& library_name)
{
  py::module importlib = py::module::import("importlib");
  return (!importlib.attr("find_loader")(library_name.c_str()).is_none());
}

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

  if (import_exists("_vxl_contrib"))
    m.attr("contrib") = py::module::import("_vxl_contrib");

}
