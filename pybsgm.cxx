#include "pybsgm.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <ios>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <bsgm/bsgm_disparity_estimator.h>
#include <bsgm/bsgm_multiscale_disparity_estimator.h>

namespace py = pybind11;

namespace pyvxl { namespace bsgm {

void wrap_bsgm(py::module &m)
{
  py::class_<bsgm_disparity_estimator_params>(m, "disparity_estimator_params")
    .def(py::init<>())
    .def_readwrite("use_16_directions", &bsgm_disparity_estimator_params::use_16_directions)
    .def_readwrite("p1_scale", &bsgm_disparity_estimator_params::p1_scale)
    .def_readwrite("p2_scale", &bsgm_disparity_estimator_params::p2_scale)
    .def_readwrite("use_gradient_weighted_smoothing", &bsgm_disparity_estimator_params::use_gradient_weighted_smoothing)
    .def_readwrite("max_grad", &bsgm_disparity_estimator_params::max_grad)
    .def_readwrite("perform_quadratic_interp", &bsgm_disparity_estimator_params::perform_quadratic_interp)
    .def_readwrite("error_check_mode", &bsgm_disparity_estimator_params::error_check_mode)
    .def_readwrite("shadow_thresh", &bsgm_disparity_estimator_params::shadow_thresh)
    .def_readwrite("census_weight", &bsgm_disparity_estimator_params::census_weight)
    .def_readwrite("xgrad_weight", &bsgm_disparity_estimator_params::xgrad_weight)
    .def_readwrite("census_tol", &bsgm_disparity_estimator_params::census_tol)
    .def_readwrite("census_rad", &bsgm_disparity_estimator_params::census_rad)
    .def_readwrite("print_timing", &bsgm_disparity_estimator_params::print_timing);

  py::class_<bsgm_multiscale_disparity_estimator>(m, "multiscale_disparity_estimator")

    .def(py::init<bsgm_disparity_estimator_params, int, int, int, int, int>(),
        py::arg("params"),py::arg("img_width"),py::arg("img_height"),
        py::arg("num_disparities"),py::arg("num_active_disparities"),py::arg("downscale_exponent"))
  /** multiscale_mode:
   * 0: constant min disparity for all pixels
   * 1: fixed min disparity in large block regions
   * 2: per-pixel min disparity based on coarse disparity
   **/
    .def("compute", [](bsgm_multiscale_disparity_estimator &est,
                       const vil_image_view<vxl_byte>& img_target,
                       const vil_image_view<vxl_byte>& img_ref,
                       const vil_image_view<bool>& invalid_target,
                       int min_disparity,
                       float invalid_disparity,
                       int multi_scale_mode=1)
         {
         vil_image_view<float> disp_target;
         bool status = est.compute(img_target, img_ref, invalid_target,
                                   min_disparity, invalid_disparity, multi_scale_mode,
                                   disp_target);
         if (!status) {
         throw std::runtime_error("bsgm_multiscale_disparity_estimator::compute() returned error");
         }
         return disp_target;
         },
         py::arg("img_target"), py::arg("img_ref"), py::arg("invalid_target"),
         py::arg("min_disparity"), py::arg("invalid_disparity"), py::arg("multi_scale_mode")
         );
}
}}
