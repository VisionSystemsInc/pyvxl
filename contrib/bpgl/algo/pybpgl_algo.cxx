#include "pybpgl_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <math.h>
#include <memory>
#include <sstream>
#include <vector>
#include <array>

#include <vgl/vgl_box_2d.h>
#include <vpgl/vpgl_affine_camera.h>
#include <vpgl/vpgl_perspective_camera.h>
#include <bpgl/algo/bpgl_heightmap_from_disparity.h>
#include <bpgl/algo/bpgl_rectify_image_pair.h>

namespace py = pybind11;

namespace pyvxl {
namespace bpgl {
namespace algo {

template <class CAMT>
void wrap_bpgl_rectify_image_pair(py::module &m, const char* name)
{
  py::class_<bpgl_rectify_image_pair<CAMT> > (m, name)
    .def(py::init<double, size_t, double, float, double, int>(),
         py::arg("min_disparity_z") = NAN, py::arg("n_points") = 1000,
         py::arg("upsample_scale") = 1.0, py::arg("invalid_pixel_val") = 0.0f,
         py::arg("min_overlap_fraction") = 0.25, py::arg("window_padding") = 0)
    /* .def("set_images", &bpgl_rectify_image_pair<CAMT>::set_images, */
    /*      py::arg("view0"), py::arg("view1")) */
    .def("set_cameras", &bpgl_rectify_image_pair<CAMT>::set_cameras,
         py::arg("cam0"), py::arg("cam1"))
    .def("set_homographies", &bpgl_rectify_image_pair<CAMT>::set_homographies,
         py::arg("H0"), py::arg("H1"), py::arg("out_ni"), py::arg("out_nj"))
    .def("compute_rectification", &bpgl_rectify_image_pair<CAMT>::compute_rectification,
         py::arg("scene_box"), py::arg("ni0") = -1, py::arg("nj0") = -1,
         py::arg("ni1") = -1, py::arg("nj1") = -1)
    .def("rectify_window", &bpgl_rectify_image_pair<CAMT>::rectify_window,
         py::arg("window"), py::arg("H"), py::arg("ni"), py::arg("nj"), py::arg("padding") = 0)
    .def("rectify_window_pair", &bpgl_rectify_image_pair<CAMT>::rectify_window_pair,
         py::arg("target_window"), py::arg("min_disparity"), py::arg("max_disparity"))
    .def("warp_image_pair", &bpgl_rectify_image_pair<CAMT>::warp_image_pair)
    .def("process", static_cast<void (bpgl_rectify_image_pair<CAMT>::*)(vgl_box_3d<double> const&, vgl_box_2d<int>&, int, int)>(&bpgl_rectify_image_pair<CAMT>::process),
         py::arg("scene_box"), py::arg("target_window") = vgl_box_2d<int>(), py::arg("min_disparity"), py::arg("max_disparity"))
    .def("process", static_cast<void (bpgl_rectify_image_pair<CAMT>::*)(vil_image_view_base_sptr const&, vil_image_view_base_sptr const&, CAMT&, CAMT&, vgl_box_3d<double>const&, vgl_box_2d<int>&, int, int)>(&bpgl_rectify_image_pair<CAMT>::process),
         py::arg("view_sptr0"), py::arg("view_sptr1"), py::arg("cam0"), py::arg("cam1"), py::arg("scene_box"), py::arg("target_window") = vgl_box_2d<int>(), py::arg("min_disparity"), py::arg("max_disparity"))
    .def("process", static_cast<void (bpgl_rectify_image_pair<CAMT>::*)(vil_image_view<unsigned char> const&, vil_image_view<unsigned char> const&, CAMT&, CAMT&, vgl_box_3d<double>const&, vgl_box_2d<int>&, int, int)>(&bpgl_rectify_image_pair<CAMT>::process),
         py::arg("view0"), py::arg("view1"), py::arg("cam0"), py::arg("cam1"), py::arg("scene_box"), py::arg("target_window") = vgl_box_2d<int>(), py::arg("min_disparity"), py::arg("max_disparity"))
    .def("process", static_cast<void (bpgl_rectify_image_pair<CAMT>::*)(vil_image_resource_sptr const&, vil_image_resource_sptr const&, CAMT&, CAMT&, vgl_box_3d<double>const&, vgl_box_2d<int>&, int, int)>(&bpgl_rectify_image_pair<CAMT>::process),
         py::arg("view0"), py::arg("view1"), py::arg("cam0"), py::arg("cam1"), py::arg("scene_box"), py::arg("target_window") = vgl_box_2d<int>(), py::arg("min_disparity"), py::arg("max_disparity"))

    .def_property_readonly("rect_window0", &bpgl_rectify_image_pair<CAMT>::rect_window0)
    .def_property_readonly("rect_window1", &bpgl_rectify_image_pair<CAMT>::rect_window1)
    .def_property_readonly("cam0", &bpgl_rectify_image_pair<CAMT>::cam0)
    .def_property_readonly("cam1", &bpgl_rectify_image_pair<CAMT>::cam1)
    /* .def_property_readonly("rect_cam0", &bpgl_rectify_image_pair<CAMT>::rect_cam0) */
    /* .def_property_readonly("rect_cam1", &bpgl_rectify_image_pair<CAMT>::rect_cam1) */
    .def_property_readonly("rect_view0", &bpgl_rectify_image_pair<CAMT>::rectified_fview0)
    .def_property_readonly("rect_view1", &bpgl_rectify_image_pair<CAMT>::rectified_fview1)
    .def_property_readonly("H0", &bpgl_rectify_image_pair<CAMT>::H0)
    .def_property_readonly("H1", &bpgl_rectify_image_pair<CAMT>::H1)
    .def_property_readonly("rectified_dims", &bpgl_rectify_image_pair<CAMT>::rectified_dims)
  ;
}

void wrap_bpgl_algo(py::module &m)
{
  m
    .def("heightmap_from_disparity", &bpgl_heightmap_from_disparity<float,  vpgl_affine_camera<double> >)
    .def("heightmap_from_disparity", &bpgl_heightmap_from_disparity<double, vpgl_affine_camera<double> >)
  ;
  wrap_bpgl_rectify_image_pair<vpgl_affine_camera<double> >(m, "rectify_image_pair_affine");
  wrap_bpgl_rectify_image_pair<vpgl_perspective_camera<double> >(m, "rectify_image_pair_perspective");
}


}}}

PYBIND11_MODULE(_bpgl_algo, m)
{
  m.doc() = "Python bindings for the VXL computer vision libraries";

  pyvxl::bpgl::algo::wrap_bpgl_algo(m);
}

