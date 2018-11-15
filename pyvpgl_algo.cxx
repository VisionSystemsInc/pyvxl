#include "pyvpgl_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vpgl/algo/vpgl_backproject.h>
#include <vpgl/algo/vpgl_camera_convert.h>
#include <vpgl/vpgl_rational_camera.h>
#include <vpgl/vpgl_affine_camera.h>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_box_3d.h>

namespace py = pybind11;

namespace pyvxl { namespace vpgl { namespace algo {

void wrap_vpgl_algo(py::module &m)
{
  py::module bproj_mod = m.def_submodule("backproject");
  bproj_mod.def("bproj_plane", [](vpgl_rational_camera<double> const& rcam,
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

  py::module aff_conv_mod = m.def_submodule("affine_camera_convert");
  aff_conv_mod.def("convert", [](vpgl_local_rational_camera<double> const& rcam,
                        vgl_box_3d<double> const& roi,
                        unsigned int num_points=10000)
  {
    vpgl_affine_camera<double> cam_out;
    bool status = vpgl_affine_camera_convert::convert(rcam, roi, cam_out, num_points);
    if (!status) {
      throw std::runtime_error("vpgl_affine_camera_convert::convert() returned error");
    }
    return cam_out;
  });

}
}}}
