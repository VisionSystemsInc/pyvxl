#include "pyvpgl_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vpgl/algo/vpgl_backproject.h>
#include <vpgl/algo/vpgl_camera_convert.h>
#include <vpgl/algo/vpgl_camera_compute.h>
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
                        unsigned int num_points)
  {
    vpgl_affine_camera<double> cam_out;
    bool status = vpgl_affine_camera_convert::convert(rcam, roi, cam_out, num_points);
    if (!status) {
      throw std::runtime_error("vpgl_affine_camera_convert::convert() returned error");
    }
    return cam_out;
  },
  "convert vpgl_local_rational_camera to vpgl_affine_camera",
  py::arg("rcam"), py::arg("roi"), py::arg("num_points")=10000);

  py::module persp_compute_mod = m.def_submodule("perspective_camera_compute");
  persp_compute_mod.def("compute_dlt", [](std::vector<vgl_point_2d<double> > const& image_pts,
                                          std::vector<vgl_point_3d<double> > const& world_pts)
  {
    vpgl_perspective_camera<double> camera;
    double err;
    bool result = vpgl_perspective_camera_compute::compute_dlt(image_pts, world_pts, camera, err);
    if(!result) {
    throw std::runtime_error("error computing perspective camera");
    }
    return camera;
  },
  "compute vpgl_perspective_camera from 2D->3D correspondences",
  py::arg("image points"), py::arg("world points"));


  py::module affine_compute_mod = m.def_submodule("affine_camera_compute");
  affine_compute_mod.def("compute", [](std::vector<vgl_point_2d<double> > const& image_pts,
                                       std::vector<vgl_point_3d<double> > const& world_pts)
  {
    vpgl_affine_camera<double> camera;
    bool result = vpgl_affine_camera_compute::compute(image_pts, world_pts, camera);
    if(!result) {
    throw std::runtime_error("error computing affine camera");
    }
    return camera;
  },
  "compute vpgl_affine_camera from 2D->3D correspondences",
  py::arg("image points"), py::arg("world points"));

}
}}}

PYBIND11_MODULE(_vpgl_algo, m)
{
  m.doc() =  "Python bindings for the VPGL Algo computer vision libraries";

  pyvxl::vpgl::algo::wrap_vpgl_algo(m);
}
