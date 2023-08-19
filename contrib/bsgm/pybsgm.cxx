#include "pybsgm.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <bsgm/bsgm_disparity_estimator.h>
#include <bsgm/bsgm_prob_pairwise_dsm.h>

#include "../../pyvxl_util.h"

namespace py = pybind11;

namespace pyvxl { namespace bsgm {

// wrapping for bsgm_prob_pairwise_dsm
// templated over both camera type and image pixel type
template<class CAM_T, class PIX_T>
void wrap_bsgm_prob_pairwise_dsm(py::module &m, std::string const& class_name)
{
  using BSGM_T = bsgm_prob_pairwise_dsm<CAM_T, PIX_T>;
  using IMAGE_T = vil_image_view<PIX_T>;

  py::class_<BSGM_T> (m, class_name.c_str())
    .def(py::init<>())
    .def(py::init<IMAGE_T const&, CAM_T const&,
         IMAGE_T const&, CAM_T const&,
         vgl_vector_3d<float> const& , vgl_vector_3d<float> const& >(),
         py::arg("view_0"), py::arg("cam_0"),
         py::arg("view_1"), py::arg("cam_1"),
         py::arg("sun_dir_0")=vgl_vector_3d<float>(0,0,0),
         py::arg("sun_dir_1")=vgl_vector_3d<float>(0,0,0)
         )
    .def("set_images_and_cams",
        overload_cast_<IMAGE_T const&, CAM_T const&,
         IMAGE_T const&, CAM_T const&, vgl_vector_3d<float> const&, vgl_vector_3d<float> const&>()
        (&BSGM_T::set_images_and_cams),
         "set images and cameras",
         py::arg("view_0"), py::arg("cam_0"),
         py::arg("view_1"), py::arg("cam_1"),
         py::arg("sun_dir_0")=vgl_vector_3d<float>(0,0,0),
         py::arg("sun_dir_1")=vgl_vector_3d<float>(0,0,0)
        )
    .def("set_dynamic_range_table", &BSGM_T::set_dynamic_range_table, py::arg("bits_per_pix_factors"),
         "the amount to scale appearance quantities with respect to effective bits per pixel")

    .def_property("params",
        overload_cast_<>()(&BSGM_T::params, py::const_),
        overload_cast_<pairwise_params const&>()(&BSGM_T::params),
        "pairwise parameters")

    .def_property("min_disparity",
        overload_cast_<>()(&BSGM_T::min_disparity, py::const_),
        overload_cast_<int>()(&BSGM_T::min_disparity),
        "minimum disparity to start search along an epipolar line")

    .def_property("max_disparity",
        overload_cast_<>()(&BSGM_T::max_disparity, py::const_),
        overload_cast_<int>()(&BSGM_T::max_disparity),
        "maximum disparity to end search along an epipolar line")

    .def_property_readonly("num_disparities",
        &BSGM_T::num_disparities,
        "number of disparities")

    .def_property_readonly("num_active_disparities",
        &BSGM_T::num_active_disparities,
        "how many disparity values are searched around the coarse search result")

    .def_property("max_height",
        overload_cast_<>()(&BSGM_T::max_height, py::const_),
        overload_cast_<float>()(&BSGM_T::max_height),
        "maximum height to define maximum and minimum disparity")

    .def_property("midpoint_z",
        overload_cast_<>()(&BSGM_T::midpoint_z, py::const_),
        overload_cast_<double>()(&BSGM_T::midpoint_z),
        "plane elevation for minimum least squares disparity")

    .def_property("scene_box",
        overload_cast_<>()(&BSGM_T::scene_box, py::const_),
        overload_cast_<vgl_box_3d<double> >()(&BSGM_T::scene_box),
        "scene box for analysis")

    .def_property("target_window",
        overload_cast_<>()(&BSGM_T::target_window, py::const_),
        overload_cast_<vgl_box_2d<int> >()(&BSGM_T::target_window),
        "Crop processing window within target image")

    /* .def_property("reference_window", */
    /*     overload_cast_<>()(&BSGM_T::reference_window, py::const_), */
    /*     overload_cast_<vgl_box_2d<int> >()(&BSGM_T::reference_window), */
    /*     "Crop processing window within reference image") */

    .def_property("H0",
        overload_cast_<>()(&BSGM_T::H0, py::const_),
        overload_cast_<vnl_matrix_fixed<double, 3, 3> >()(&BSGM_T::H0),
        "Image 0's rectification homography")

    .def_property("H1",
        overload_cast_<>()(&BSGM_T::H1, py::const_),
        overload_cast_<vnl_matrix_fixed<double, 3, 3> >()(&BSGM_T::H1),
        "Image 1's rectification homography")

    .def_property("rect_ni",
        overload_cast_<>()(&BSGM_T::rect_ni, py::const_),
        overload_cast_<size_t>()(&BSGM_T::rect_ni),
        "Width of rectified images")

    .def_property("rect_nj",
        overload_cast_<>()(&BSGM_T::rect_nj, py::const_),
        overload_cast_<size_t>()(&BSGM_T::rect_nj),
        "Height of rectified images")

    .def("rectified_bview0", &BSGM_T::rectified_bview0,
         "rectified image view 0")
    .def("rectified_bview1", &BSGM_T::rectified_bview1,
         "rectified image view 1")

    .def("rectified_cam0", &BSGM_T::rectified_cam0,
         "rectified camera 0")
    .def("rectified_cam1", &BSGM_T::rectified_cam1,
         "rectified camera 1")

    .def("invalid_map_fwd", &BSGM_T::invalid_map_fwd,
         "invalid map for forward disparity")
    .def("invalid_map_rev", &BSGM_T::invalid_map_rev,
         "invalid map for reverse disparity")

    .def("disparity_fwd", &BSGM_T::disparity_fwd,
         "forward disparity")
    .def("disparity_rev", &BSGM_T::disparity_rev,
         "reverse disparity")

    .def("tri_3d_fwd", &BSGM_T::tri_3d_fwd,
         "forward triangulation result")
    .def("tri_3d_rev", &BSGM_T::tri_3d_rev,
         "reverse triangulation result")
    .def("xyz_prob", &BSGM_T::xyz_prob,
         "image with xyz and probability")

    .def("heightmap_fwd", &BSGM_T::heightmap_fwd,
         "forward heightmap")
    .def("heightmap_rev", &BSGM_T::heightmap_rev,
         "reverse heightmap")

    .def("ptset_fwd", &BSGM_T::ptset_fwd,
         "forward pointset")
    .def("ptset_rev", &BSGM_T::ptset_rev,
         "reverse pointset")

    .def("prob_ptset", &BSGM_T::prob_ptset,
         "probabilistic pointset")
    .def("prob_pdf", &BSGM_T::prob_pdf,
         "probabilistic pdf")

    .def("prob_heightmap", &BSGM_T::prob_heightmap,
         "probabilistic heightmap")
    .def("prob_confidence", &BSGM_T::prob_confidence,
         "probabilistic confidence")
    .def("radial_std_dev_image", &BSGM_T::radial_std_dev_image,
         "radial standard deviation")

    .def("rect_target_stype", &BSGM_T::rect_target_stype,
         "surface types in rectified target space")
    .def("dsm_grid_stype", &BSGM_T::dsm_grid_stype,
         "surface types in dsm space")

    .def("rectify", &BSGM_T::rectify,
         py::call_guard<py::gil_scoped_release>(),
         "image rectification")

    .def("compute_disparity_fwd", &BSGM_T::compute_disparity_fwd,
         py::call_guard<py::gil_scoped_release>(),
         "compute forward disparity")
    .def("compute_disparity_rev", &BSGM_T::compute_disparity_rev,
         py::call_guard<py::gil_scoped_release>(),
         "compute reverse disparity")

    .def("compute_height_fwd", &BSGM_T::compute_height_fwd,
         py::call_guard<py::gil_scoped_release>(),
         "compute forward height data (tri_3d, ptset, heightmap)")
    .def("compute_height_rev", &BSGM_T::compute_height_rev,
         py::call_guard<py::gil_scoped_release>(),
         "compute reverse height data (tri_3d, ptset, heightmap)")

    .def("compute_prob", &BSGM_T::compute_prob,
         py::call_guard<py::gil_scoped_release>(),
         "compute probabilistic height")

    .def("process", &BSGM_T::process,
         py::call_guard<py::gil_scoped_release>(),
         py::arg("with_consistency_check") = true,
         py::arg("knn_consistency") = true,
         py::arg("compute_fwd_rev_ptsets_hmaps") = true,
         "Main process method")

    .def("process_with_windows", &BSGM_T::process_with_windows,
         py::call_guard<py::gil_scoped_release>(),
         py::arg("first_window") = true,
         py::arg("with_consistency_check")= false,
         py::arg("print_timing") = false,
         "Main process method for windows into stereo image pair")

    .def("save_prob_ptset_color", &BSGM_T::save_prob_ptset_color,
         py::arg("path"),
         "apply a color map to probabilty values and output an ascii color point cloud")

    .def("save_rect_target_stype", &BSGM_T::save_rect_target_stype,
         py::arg("path"),
         "save surface types in rectified target space")

    .def("save_dsm_grid_stype", &BSGM_T::save_dsm_grid_stype,
         py::arg("path"),
         "save surface types in dsm space")

    ;
}



void wrap_bsgm(py::module &m)
{

  // bsgm_disparity_estimator_params
  py::class_<bsgm_disparity_estimator_params> (m, "disparity_estimator_params")
    .def(py::init(&init_struct_from_kwargs<bsgm_disparity_estimator_params>))
    .def("__repr__", repr_by_dict<bsgm_disparity_estimator_params>)
    .def("as_dict", struct_to_dict<bsgm_disparity_estimator_params>)
    .def("set", set_struct_from_kwargs<bsgm_disparity_estimator_params>)
    .def_readwrite("use_16_directions", &bsgm_disparity_estimator_params::use_16_directions,
                   "Use 16 directions in the dynamic programming, otherwise 8")
    .def_readwrite("p1_scale", &bsgm_disparity_estimator_params::p1_scale,
                   "Scale the internally set P1 smoothing parameter")
    .def_readwrite("p2_scale", &bsgm_disparity_estimator_params::p2_scale,
                   "Scale the internally set P2 smoothing parameter")
    .def_readwrite("use_gradient_weighted_smoothing", &bsgm_disparity_estimator_params::use_gradient_weighted_smoothing,
                   "Use gradient-weighted P2 smoothing")
    .def_readwrite("use_shadow_step_p2_adjustment", &bsgm_disparity_estimator_params::use_shadow_step_p2_adjustment,
                   "Adjust P2 according to probability of shadow step")
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
    .def(py::init(&init_struct_from_kwargs<pairwise_params>))
    .def("__repr__", repr_by_dict<pairwise_params>)
    .def("as_dict", struct_to_dict<pairwise_params>)
    .def("set", set_struct_from_kwargs<pairwise_params>)
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
    .def_readwrite("min_overlap_fraction", &pairwise_params::min_overlap_fraction_,
                   "minimum fraction of points in the zero disparity plane (w/in scene box) that must project into both images")
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
    .def_readwrite("coarse_dsm_disparity_estimate", &pairwise_params::coarse_dsm_disparity_estimate_,
                   "use the reduced resolution dsm to estimate min disparity")
    .def_readwrite("effective_bits_per_pixel", &pairwise_params::effective_bits_per_pixel_,
                   "The actual intensity dynamic range, e.g. 11 bits")
    .def_readwrite("window_padding", &pairwise_params::window_padding_,
                   "How many pixels to pad the target window by, post rectification")
    ;

  // bsgm_prob_pairwise_dsm
  wrap_bsgm_prob_pairwise_dsm<vpgl_affine_camera<double>, unsigned char >(m, "prob_pairwise_dsm_affine");//legacy default is byte
  wrap_bsgm_prob_pairwise_dsm<vpgl_affine_camera<double>, unsigned short >(m, "prob_pairwise_dsm_affine_short");
  wrap_bsgm_prob_pairwise_dsm<vpgl_perspective_camera<double>, unsigned char>(m, "prob_pairwise_dsm_perspective");//legacy default is byte
  wrap_bsgm_prob_pairwise_dsm<vpgl_perspective_camera<double>, unsigned short>(m, "prob_pairwise_dsm_perspective_short");

} // wrap_bsgm

}}


PYBIND11_MODULE(_bsgm, m)
{
  m.doc() =  "Python bindings for the VXL BSGM computer vision library";
  pyvxl::bsgm::wrap_bsgm(m);
}
