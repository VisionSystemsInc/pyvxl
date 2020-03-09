#include "pybsgm.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <bsgm/bsgm_disparity_estimator.h>
#include <bsgm/bsgm_prob_pairwise_dsm.h>

namespace py = pybind11;

namespace pyvxl { namespace bsgm {

void wrap_bsgm(py::module &m)
{

  // bsgm_disparity_estimator_params
  py::class_<bsgm_disparity_estimator_params> (m, "disparity_estimator_params")
    .def(py::init<>())
    .def_readwrite("use_16_directions", &bsgm_disparity_estimator_params::use_16_directions,
                   "Use 16 directions in the dynamic programming, otherwise 8")
    .def_readwrite("p1_scale", &bsgm_disparity_estimator_params::p1_scale,
                   "Scale the internally set P1 smoothing parameter")
    .def_readwrite("p2_scale", &bsgm_disparity_estimator_params::p2_scale,
                   "Scale the internally set P2 smoothing parameter")
    .def_readwrite("use_gradient_weighted_smoothing", &bsgm_disparity_estimator_params::use_gradient_weighted_smoothing,
                   "Use gradient-weighted P2 smoothing")
    .def_readwrite("max_grad", &bsgm_disparity_estimator_params::max_grad,
                   "In gradient-weighted smoothing, gradients beyond this magnitude are truncated")
    .def_readwrite("perform_quadratic_interp", &bsgm_disparity_estimator_params::perform_quadratic_interp,
                   "Use quadratic interpolation to obtain sub-pixel estimates of final disparity map")
    .def_readwrite("error_check_mode", &bsgm_disparity_estimator_params::error_check_mode,
                   "Mode for finding and fixing errors in the disparity map (0 = raw disparity map, 1 = flag bad pixels, 2 = interpolate bad pixels)")
    .def_readwrite("shadow_thresh", &bsgm_disparity_estimator_params::shadow_thresh,
                   "When set > 0, pixels below this threshold will be flagged as invalid when error_check_mode > 0")
    .def_readwrite("bias_weight", &bsgm_disparity_estimator_params::bias_weight,
                   "Strength of SGM directional average bias")
    .def_readwrite("bias_dir", &bsgm_disparity_estimator_params::bias_dir,
                   "Direction of SGM directional average bias")
    .def_readwrite("census_weight", &bsgm_disparity_estimator_params::census_weight,
                   "Census appearance cost weighting")
    .def_readwrite("xgrad_weight", &bsgm_disparity_estimator_params::xgrad_weight,
                   "X-gradient appearance cost weighting")
    .def_readwrite("census_tol", &bsgm_disparity_estimator_params::census_tol,
                   "Pixel differences less than this magnitude are not considered in the census computation")
    .def_readwrite("census_rad", &bsgm_disparity_estimator_params::census_rad,
                   "Length (1,2, or 3) of the census kernal will be 2*census_rad+1")
    .def_readwrite("print_timing", &bsgm_disparity_estimator_params::print_timing,
                   "Print detailed timing information to cerr")
    ;

  // bsgm_prob_pairwise_dsm pairwise_params
  py::class_<pairwise_params> (m, "pairwise_params")
    .def(py::init<>())
    .def_readwrite("de_params", &pairwise_params::de_params_,
                   "disparity estimator params")
    .def_property("shadow_thresh",
                  [](pairwise_params& self) { return self.shadow_thresh_; },
                  &pairwise_params::shadow_thresh,
                  "intensity level out of 255 below which is considered to be in shadow, thus invalid")
    .def_property("quad_interp",
                  [](pairwise_params& self) { return self.quad_interp_; },
                  &pairwise_params::quad_interp,
                  "if true, perform quadratic interpolation of disparity with respect to cost")
    .def_readwrite("multi_scale_mode", &pairwise_params::multi_scale_mode_,
                   "disparity estimator multi-scale mode")
    .def_readwrite("active_disparity_factor", &pairwise_params::active_disparity_factor_,
                   "fraction of full disparity range is used for fine search")
    .def_readwrite("downscale_exponent", &pairwise_params::downscale_exponent_,
                   "in coarse to fine disparity, what is the downsample ratio as 2^exponent")
    .def_readwrite("ground_sample_dist", &pairwise_params::ground_sample_dist_,
                   "height map grid spacing, also relates to consistent distance tolerance")
    .def_readwrite("upsample_scale_factor", &pairwise_params::upsample_scale_factor_,
                   "upsample the rectified images by scale factor")
    .def_readwrite("std_dev", &pairwise_params::std_dev_,
                   "standard deviation of consistent disparity point distances")
    .def_readwrite("use_z_vs_d_prob", &pairwise_params::use_z_vs_d_prob_,
                   "multiply height probabilty with additional z vs d scale probability factor")
    .def_readwrite("min_z_vs_d_scale", &pairwise_params::min_z_vs_d_scale_,
                   "lowest z vs d scale factor that is typically obtained in meters/pixel")
    .def_readwrite("z_vs_d_std_dev", &pairwise_params::z_vs_d_std_dev_,
                   "standard deviation for the z vs d Gaussian distribution")
    .def_readwrite("min_neighbors", &pairwise_params::min_neighbors_,
                   "pointset->heightmap gridding minimum number of neighbors")
    .def_readwrite("max_neighbors", &pairwise_params::max_neighbors_,
                   "pointset->heightmap gridding maximum number of neighbors")
    .def_readwrite("neighbor_dist_factor", &pairwise_params::neighbor_dist_factor_,
                   "pointset->heightmap gridding distance of (neighbor_dist_factor * ground_sample_distance)")
    ;

} // wrap_bsgm

}}


PYBIND11_MODULE(_bsgm, m)
{
  m.doc() =  "Python bindings for the VXL BSGM computer vision library";
  pyvxl::bsgm::wrap_bsgm(m);
}
