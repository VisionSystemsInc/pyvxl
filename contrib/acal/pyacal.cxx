#include "pyacal.h"
#include "../../pyvxl_util.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <acal/acal_f_utils.h>
#include <acal/acal_match_graph.h>
#include <acal/acal_match_utils.h>

// io classes for py::pickle
#include <acal/io/acal_io_match_graph.h>
#include <acal/io/acal_io_match_utils.h>

namespace py = pybind11;

namespace pyvxl { namespace acal {

void wrap_acal(py::module &m)
{

  // acal_f_utils::f_params
  py::class_<f_params> (m, "f_params")
    .def(py::init(&init_struct_from_kwargs<f_params>))
    .def("__repr__", repr_by_dict<f_params>)
    .def("as_dict", struct_to_dict<f_params>)
    .def("set", set_struct_from_kwargs<f_params>)
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

  // acal_match_utils::acal_corr
  py::class_<acal_corr>(m, "corr")
    .def(py::init(&init_struct_from_kwargs<acal_corr>))
    .def("__repr__", repr_by_dict<acal_corr>)
    .def("as_dict", struct_to_dict<acal_corr>)
    .def("set", set_struct_from_kwargs<acal_corr>)
    .def_readwrite("id", &acal_corr::id_)
    .def_readwrite("pt", &acal_corr::pt_)
    .def(py::self == py::self)
    .def(py::pickle(&vslPickleGetState<acal_corr>,
                    &vslPickleSetState<acal_corr>))
    ;

  // acal_match_utils::acal_match_pair
  py::class_<acal_match_pair>(m, "match_pair")
    .def(py::init(&init_struct_from_kwargs<acal_match_pair>))
    .def("__repr__", repr_by_dict<acal_match_pair>)
    .def("as_dict", struct_to_dict<acal_match_pair>)
    .def("set", set_struct_from_kwargs<acal_match_pair>)
    .def_readwrite("corr1", &acal_match_pair::corr1_)
    .def_readwrite("corr2", &acal_match_pair::corr2_)
    .def(py::self == py::self)
    .def(py::pickle(&vslPickleGetState<acal_match_pair>,
                    &vslPickleSetState<acal_match_pair>))
    ;

  // match_params
  py::class_<match_params> (m, "match_params")
    .def(py::init(&init_struct_from_kwargs<match_params>))
    .def("__repr__", repr_by_dict<match_params>)
    .def("as_dict", struct_to_dict<match_params>)
    .def("set", set_struct_from_kwargs<match_params>)
    .def_readwrite("min_n_tracks", &match_params::min_n_tracks_,
                   "minimum number of tracks for a graph edge")
    .def_readwrite("min_n_cams", &match_params::min_n_cams_,
                   "minimum number of cameras in a graph clique")
    .def_readwrite("max_proj_error", &match_params::max_proj_error_,
                   "max projection error for a reasonable solution")
    .def_readwrite("max_uncal_proj_error", &match_params::max_uncal_proj_error_,
                   "max initial projection error")
    .def(py::self == py::self)
    .def(py::pickle(&vslPickleGetState<match_params>,
                    &vslPickleSetState<match_params>))
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

