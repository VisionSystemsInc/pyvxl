#include <pybind11/pybind11.h>
#include "pyvnl.h"
#include "pyvgl.h"
#include "pyvpgl.h"
#include "pyvil.h"
#include "pyvgl_algo.h"

namespace py = pybind11;

PYBIND11_MODULE(vxl, m)
{
  m.doc() =  "Python bindings for the VXL computer vision libraries";

  pyvxl::wrap_vnl(m);
  pyvxl::wrap_vgl(m);
  pyvxl::wrap_vpgl(m);
  pyvxl::wrap_vil(m);
  pyvxl::wrap_vgl_algo(m);
}
