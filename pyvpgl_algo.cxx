#include "pyvpgl_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vpgl/algo/vpgl_backproject.h>
#include <vpgl/vpgl_rational_camera.h>
#include <vgl/vgl_plane_3d.h>

namespace py = pybind11;

namespace pyvxl { namespace vpgl { namespace algo {

void wrap_vpgl_algo(py::module &m)
{
  py::module mod = m.def_submodule("backproject");
  mod.def("bproj_plane", [](vpgl_rational_camera<double> const& rcam,
                            vgl_point_2d<double> const& image_point,
                            vgl_plane_3d<double> const& plane,
                            vgl_point_3d<double> const& initial_guess,
                            double error_tol = 0.05,
                            double relative_diameter = 1.0)
          {
          vgl_point_3d<double> output;
          bool status = vpgl_backproject::bproj_plane(rcam, image_point, plane, initial_guess, output, error_tol, relative_diameter);
          if (!status) {
          throw std::runtime_error("vpgl_backproject::bproj_plane() returned error");
          }
          return output;
          });
}
}}}
