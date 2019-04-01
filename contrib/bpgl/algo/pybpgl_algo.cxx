#include "pybpgl_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>
#include <vector>
#include <array>

#include <vpgl/vpgl_affine_camera.h>
#include <bpgl/algo/bpgl_heightmap_from_disparity.h>

namespace py = pybind11;

namespace pyvxl { namespace bpgl { namespace algo {

void wrap_bpgl_algo(py::module &m)
{
  m.def("heightmap_from_disparity_affine", &bpgl_heightmap_from_disparity<vpgl_affine_camera<double> >);

}

}}}

PYBIND11_MODULE(_bpgl_algo, m)
{
  m.doc() =  "Python bindings for the VXL computer vision libraries";

  pyvxl::bpgl::algo::wrap_bpgl_algo(m);

/* #ifdef PYVXL_WITH_CONTRIB_BPGL */
/*   mod = m.def_submodule("bpgl"); */
/*   mod = mod.def_submodule("algo"); */
/*   pyvxl::bpgl::algo::wrap_bpgl_algo(mod); */
/* #endif */
}
