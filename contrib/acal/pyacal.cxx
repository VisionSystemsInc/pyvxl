#include "pyacal.h"
#include "../../pyvxl_util.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <acal/acal_f_utils.h>
#include <acal/acal_match_graph.h>
#include <acal/acal_match_utils.h>

// io classes for py::pickle
#include <acal/io/acal_io_f_utils.h>
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
    .def(py::self == py::self)
    .def(py::pickle(&vslPickleGetState<f_params>,
                    &vslPickleSetState<f_params>))
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

  // acal_match_graph::match_params
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

  // acal_match_tree::acal_match_node
  py::class_<acal_match_node, std::shared_ptr<acal_match_node> >(m, "match_node")
    .def(py::init<>())
    .def("__len__", &acal_match_node::size)
    .def("is_leaf", &acal_match_node::is_leaf)
    .def("is_root", &acal_match_node::is_root)
    .def_readonly("cam_id", &acal_match_node::cam_id_)
    .def_readonly("node_depth", &acal_match_node::node_depth_)
    .def_readonly("children", &acal_match_node::children_)
    .def_readonly("self_to_child_matches", &acal_match_node::self_to_child_matches_)
    .def("parent", &acal_match_node::parent, py::arg("root"))
    ;

  // acal_match_tree::acal_match_tree
  py::class_<acal_match_tree, std::shared_ptr<acal_match_tree> >(m, "match_tree")
    .def(py::init<>())
    .def("__len__", &acal_match_tree::size)
    .def_readonly("min_n_tracks", &acal_match_tree::min_n_tracks_)
    .def_readonly("root", &acal_match_tree::root_)
    .def("save_tree_dot_format", &acal_match_tree::save_tree_dot_format,
         "save a match tree to a dot file",
         py::arg("path"))
    ;

  // acal_match_graph::match_vertex
  py::class_<match_vertex, std::shared_ptr<match_vertex> >(m, "match_vertex")
    .def(py::init<>())
    .def(py::init<size_t>())
    .def_readwrite("cam_id", &match_vertex::cam_id_)
    .def_readwrite("mark", &match_vertex::mark_)
    ;

  // acal_match_graph::match_edge
  py::class_<match_edge, std::shared_ptr<match_edge> >(m, "match_edge")
    .def(py::init<>())
    .def_readwrite("id", &match_edge::id_)
    .def_readwrite("matches", &match_edge::matches_)
    .def_readwrite("v0", &match_edge::v0_)
    .def_readwrite("v1", &match_edge::v1_)
    ;

  // acal_match_graph::acal_match_graph
  py::class_<acal_match_graph>(m, "match_graph")

    // Constructors
    .def(py::init<>())
    .def(py::init<std::map<size_t, std::map<size_t, std::vector<acal_match_pair> > > const&>())

    // Properties
    .def_property("params", &acal_match_graph::get_params, &acal_match_graph::set_params)
    .def_property("image_paths", &acal_match_graph::get_image_paths, &acal_match_graph::set_image_paths)
    .def_property("acams", &acal_match_graph::all_acams, &acal_match_graph::set_all_acams)
    .def_property("vertices", &acal_match_graph::vertices, &acal_match_graph::set_vertices)
    .def_property("edges", &acal_match_graph::edges, &acal_match_graph::set_edges)
    .def_property("connected_components", &acal_match_graph::get_connected_components, &acal_match_graph::set_connected_components)
    .def_property("focus_tracks", &acal_match_graph::get_focus_tracks, &acal_match_graph::set_focus_tracks)
    .def_property("focus_track_metrics", &acal_match_graph::get_focus_track_metrics, &acal_match_graph::set_focus_track_metrics)
    .def_property("trees", &acal_match_graph::get_match_trees, &acal_match_graph::set_match_trees)
    .def_property("tree_metrics", &acal_match_graph::get_match_tree_metrics, &acal_match_graph::set_match_tree_metrics)

    // Methods
    .def("find_connected_components", &acal_match_graph::find_connected_components,
         "Construct connected components from vertices")
    .def("compute_focus_tracks", &acal_match_graph::compute_focus_tracks,
         "Identify consistent correspondence tracks")
    .def("load_affine_cams", &acal_match_graph::load_affine_cams, 
         "Load uncorrected cameras", py::arg("affine_cam_path"))
    .def("compute_match_trees", &acal_match_graph::compute_match_trees,
         "For each focus vertex, create a tree of consistent matches")
    .def("validate_match_trees_and_set_metric", &acal_match_graph::validate_match_trees_and_set_metric,
         "validate match trees and set metric")
    .def("save_graph_dot_format", &acal_match_graph::save_graph_dot_format, 
         "save a match graph to a dot file", py::arg("path"))

    ;

} // wrap_acal

}}


PYBIND11_MODULE(_acal, m)
{
  m.doc() =  "Python bindings for the VXL ACAL computer vision library";
  pyvxl::acal::wrap_acal(m);
}

