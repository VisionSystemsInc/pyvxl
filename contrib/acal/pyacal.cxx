#include "pyacal.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <bpgl/acal/acal_f_utils.h>
#include <bpgl/acal/acal_match_graph.h>

namespace py = pybind11;

namespace pyvxl { namespace acal {

void wrap_acal(py::module &m)
{

  // f_params
  py::class_<f_params> (m, "f_params")
    .def(py::init<>())
    .def_readwrite("epi_dist_mul", &f_params::epi_dist_mul_,
                   "multiplier on the error for the lowest 10 percent of epipolar line distances")
    .def_readwrite("max_epi_dist", &f_params::max_epi_dist_,
                   "an absolute threshold on epipolar line distances in case error estimation fails")
    .def_readwrite("F_similar_abcd_tol", &f_params::F_similar_abcd_tol_,
                   "max abs difference |a+c| + |b+d|, a measure of viewpoint similarity")
    .def_readwrite("F_similar_e_tol", &f_params::F_similar_e_tol_,
                   "max abs value of offset, e to determine similar images")
    .def_readwrite("ray_uncertainty_tol", &f_params::ray_uncertainty_tol_,
                   "max ray uncertainty to keep camera pair")
    .def_readwrite("min_num_matches", &f_params::min_num_matches_,
                   "minimum number of required matches to output to fmatches file")
    ;

  // match_params
  py::class_<match_params> (m, "match_params")
    .def(py::init<>())
    .def_readwrite("min_n_tracks", &match_params::min_n_tracks_,
                   "minimum number of tracks for a graph edge")
    .def_readwrite("min_n_cams", &match_params::min_n_cams_,
                   "minimum number of cameras in a graph clique")
    .def_readwrite("max_proj_error", &match_params::max_proj_error_,
                   "max projection error for a reasonable solution")
    .def_readwrite("max_uncal_proj_error", &match_params::max_uncal_proj_error_,
                   "max initial projection error")
    ;

  // acal_match_graph
  py::class_<acal_match_graph>(m, "match_graph")
    .def(py::init<>())
    ;

} // wrap_acal

}}


PYBIND11_MODULE(_acal, m)
{
  m.doc() =  "Python bindings for the VXL ACAL computer vision library";
  pyvxl::acal::wrap_acal(m);
}

