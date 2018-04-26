#include "pyvpgl_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vpgl/algo/vpgl_camera_homographies.h>

namespace py = pybind11;

namespace pyvxl { namespace vpgl { namespace algo {

void wrap_vpgl_algo(py::module &m) {
  py::class_<vpgl_camera_homographies>(m, "camera_homographies")
      .def_static("homography_from_camera",
                  (vgl_h_matrix_2d<double> (*)(vpgl_proj_camera<double> const &,
                                               vgl_plane_3d<double> const &)) &vpgl_camera_homographies::homography_from_camera)
      .def_static("homography_from_camera",
                  (vgl_h_matrix_2d<double> (*)(vpgl_perspective_camera<double> const &,
                                               vgl_plane_3d<double> const &)) &vpgl_camera_homographies::homography_from_camera)
      .def_static("homography_to_camera",
                  (vgl_h_matrix_2d<double> (*)(vpgl_proj_camera<double> const &,
                                               vgl_plane_3d<double> const &)) &vpgl_camera_homographies::homography_to_camera)
      .def_static("homography_to_camera",
                  (vgl_h_matrix_2d<double> (*)(vpgl_perspective_camera<double> const &,
                                               vgl_plane_3d<double> const &)) &vpgl_camera_homographies::homography_to_camera)
      .def_static("transform_camera_to_plane",
                  (vpgl_proj_camera<double> (*)(vpgl_proj_camera<double> const &,
                                                vgl_plane_3d<double> const &)) &vpgl_camera_homographies::transform_camera_to_plane)
      .def_static("transform_camera_to_plane",
                  (vpgl_perspective_camera<double> (*)(
                      vpgl_perspective_camera<double> const &,
                      vgl_plane_3d<double> const &)) &vpgl_camera_homographies::transform_camera_to_plane);

}
}}}
