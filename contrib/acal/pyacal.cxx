#include "pyacal.h"
#include "py_struct.h"

#include <map>
#include <sstream>
#include <string>

#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <bpgl/acal/acal_f_utils.h>
#include <bpgl/acal/acal_match_graph.h>
#include <vpgl/vpgl_affine_camera.h>

namespace py = pybind11;

using namespace py::literals;

namespace pyvxl { namespace acal {


void wrap_f_params(py::module &m) {
  py::class_<f_params> (m, "f_params")
    .def(py::init<>())
    .def(py::init(&dict_to_struct<f_params>))
    .def("as_dict", struct_to_dict<f_params>)
    .def("__repr__", repr_by_dict<f_params>)
    .def_readwrite("epi_dist_mul", &f_params::epi_dist_mul_,
                   "multiplier on the error for the lowest 10% of epipolar line distances")
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
}


void wrap_match_params(py::module &m) {
  py::class_<match_params> (m, "match_params")
    .def(py::init<>())
    .def(py::init(&dict_to_struct<match_params>))
    .def("as_dict", struct_to_dict<match_params>)
    .def("__repr__", repr_by_dict<match_params>)
    .def_readwrite("min_n_tracks", &match_params::min_n_tracks_,
                   "minimum number of tracks for a graph edge")
    .def_readwrite("min_n_cams", &match_params::min_n_cams_,
                   "minimum number of cameras in a graph clique")
    .def_readwrite("max_proj_error", &match_params::max_proj_error_,
                   "max projection error for a reasonable solution")
    .def_readwrite("max_uncal_proj_error", &match_params::max_uncal_proj_error_,
                   "max initial projection error")
    .def(py::pickle(
        [](const match_params &mp) {  // __getstate__
          /* Return a tuple that fully encodes the state of the object */
          return py::make_tuple(mp.min_n_tracks_, mp.min_n_cams_,
                                mp.max_proj_error_, mp.max_uncal_proj_error_);
        },

        [](py::tuple t) {  // __setstate__
          if (t.size() != 4)
            throw std::runtime_error("Can't unpickle match_params - Invalid state!");

          /* Create a new C++ instance */
          match_params mp;

          /* Assign additional state */
          mp.min_n_tracks_ = t[0].cast<size_t>();
          mp.min_n_cams_ = t[1].cast<size_t>();
          mp.max_proj_error_ = t[2].cast<double>();
          mp.max_uncal_proj_error_ = t[3].cast<double>();

          return mp;
        }))
    ;
}


void wrap_acal_corr(py::module &m) {
  py::class_<acal_corr>(m, "corr")
    .def(py::init<>())
    .def(py::init<size_t, vgl_point_2d<double> const&>())
    .def_readwrite("id", &acal_corr::id_)
    .def_readwrite("pt", &acal_corr::pt_)
    .def(py::pickle(
        [](acal_corr &ac) {  // __getstate__
          // Return a tuple that fully encodes the state of the object
          return py::make_tuple(ac.id_, ac.pt_);
        },
        [](py::tuple t) {  // __setstate__
          if (t.size() != 2)
            throw std::runtime_error("Can't unpickle acal_corr - Invalid state!");

          /* Create a new C++ instance */
          acal_corr ac(t[0].cast<size_t>(), t[1].cast<vgl_point_2d<double> >());

          return ac;
        }))
    ;
}


void wrap_acal_match_pair(py::module &m) {
  py::class_<acal_match_pair>(m, "match_pair")
    .def(py::init<>())
    .def(py::init<acal_corr const&, acal_corr const&>())
    .def_readwrite("corr1", &acal_match_pair::corr1_)
    .def_readwrite("corr2", &acal_match_pair::corr2_)
    .def(py::pickle(
        [](acal_match_pair &amp) {  // __getstate__
          // Return a tuple that fully encodes the state of the object
          return py::make_tuple(amp.corr1_, amp.corr2_);
        },
        [](py::tuple t) {  // __setstate__
          if (t.size() != 2)
            throw std::runtime_error("Can't unpickle acal_match_pair - Invalid state!");

          /* Create a new C++ instance */
          acal_match_pair amp(t[0].cast<acal_corr>(), t[1].cast<acal_corr>());

          return amp;
        }))
    ;
}


void wrap_match_vertex(py::module &m) {
  py::class_<match_vertex, std::shared_ptr<match_vertex> /* <- holder type */>(m, "match_vertex")
    .def(py::init<>())
    .def(py::init<size_t>())
    .def_readwrite("cam_id", &match_vertex::cam_id_)
    .def_readwrite("mark", &match_vertex::mark_)
    ;
}


void wrap_match_edge(py::module &m) {
  py::class_<match_edge, std::shared_ptr<match_edge> /* <- holder type */>(m, "match_edge")
    .def(py::init<>())
    .def_readwrite("id", &match_edge::id_)
    .def_readwrite("matches", &match_edge::matches_)
    .def_readwrite("v0", &match_edge::v0_)
    .def_readwrite("v1", &match_edge::v1_)
    ;
}


void wrap_acal_match_tree(py::module &m) {
  py::class_<acal_match_tree, std::shared_ptr<acal_match_tree> /* <- holder type */>(m, "match_tree")
    .def(py::init<>())
    .def("save_tree_dot_format", &acal_match_tree::save_tree_dot_format, "save a match tree to a dot file",
          py::arg("path"))
    ;
}


void wrap_acal_match_graph(py::module &m) {
  py::class_<acal_match_graph>(m, "match_graph")

    // Constructors
    .def(py::init<>())

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
    .def("save_graph_dot_format", &acal_match_graph::save_graph_dot_format, "save a match graph to a dot file",
          py::arg("path"))

    // Pickle
    .def(py::pickle(
        [](acal_match_graph &mg) {  // __getstate__

          // Get graph attributes
          match_params params = mg.get_params();
          std::map<size_t, std::string>& image_paths = mg.get_image_paths();
          std::map<size_t, vpgl_affine_camera<double> > all_acams = mg.all_acams();
          std::map<size_t, std::shared_ptr<match_vertex> > vertices = mg.vertices();
          std::vector<std::shared_ptr<match_edge> > edges = mg.edges();
          std::vector<std::vector<std::shared_ptr<match_vertex> > > connected_components = mg.get_connected_components();
          std::map<size_t, std::map<size_t, std::vector< std::map<size_t, vgl_point_2d<double> > > > > focus_tracks = mg.get_focus_tracks();
          std::vector<double> focus_track_metrics = mg.get_focus_track_metrics();
          std::map<size_t, std::map<size_t, std::shared_ptr<acal_match_tree> > > match_trees = mg.get_match_trees();
          std::vector<size_t> match_tree_metrics = mg.get_match_tree_metrics();

          /* We can't (easily) directly serialize these shared pointers.
           * So instead we'll assign an ID to each object, and replace
           * references to these objects with their IDs. Then when we deserialize,
           * we'll rebuild the lower level objects and fill in the references.
           */

          // Create serializable representations of each vertex
          std::map<size_t, std::pair<bool, std::vector<size_t> > > serializable_vertices;
          for (auto const& item : vertices) {
            size_t vertex_id = item.first;
            std::shared_ptr<match_vertex> vertex = item.second;

            size_t cam_id = vertex->cam_id_;
            if (serializable_vertices.count(cam_id) > 0) {
              throw std::runtime_error("Can't pickle acal_match_graph - Non-unique vertex cam ID!");
            }

            // vector of edges to vector of edge IDs
            std::vector<size_t> edge_ids;
            for (match_edge* me_ptr : vertex->edges_) {
              edge_ids.push_back(me_ptr->id_);
            }

            serializable_vertices[cam_id] = std::make_pair(vertex->mark_,
                                                            edge_ids);
          }

          // Create serializable representations of each edge
          std::map<size_t, std::tuple<size_t, size_t, std::vector<acal_match_pair> > > serializable_edges;
          for (auto const& e_ptr : edges) {
            size_t edge_id = e_ptr->id_;
            if (serializable_edges.count(edge_id) > 0) {
              throw std::runtime_error("Can't pickle acal_match_graph - Non-unique edge ID!");
            }
            serializable_edges[edge_id] = std::make_tuple(e_ptr->v0_->cam_id_,
                                                          e_ptr->v1_->cam_id_,
                                                          e_ptr->matches_);
          }

          // Replace the connected component vertices with IDs
          std::vector<std::vector<size_t> > connected_components_with_ids;
          for (auto const& cc : connected_components) {
            std::vector<size_t> verts_in_cc;
            for (auto const& v_ptr : cc) {
              verts_in_cc.push_back(v_ptr->cam_id_);
            }
            connected_components_with_ids.push_back(verts_in_cc);
          }

          // Collect all the acal_match_node objects across all the trees
          // I'm doing this because I'm not sure if the different trees share nodes
          std::vector<acal_match_node*> all_nodes;
          for (auto const& item1 : match_trees) {
            for (auto const& item2 : item1.second) {
              auto tree = item2.second;
              tree->collect_nodes(tree->root_, all_nodes);
            }
          }

          // Remove duplicates (Actually, I don't think this is necessary since the below will just
          // overwrite a duplicate object with a new, still unique, id)

          // Create an "inverse" mapping, from acal_match_tree objects in memory, to new node IDs
          std::map<acal_match_node*, int> node_id_map;
          for (int node_id=0; node_id < all_nodes.size(); node_id++) {
            node_id_map[all_nodes[node_id]] = node_id;
          }

          // Use these new node IDs to create serializable representations of each node
          std::map<int, std::tuple<size_t, size_t, int, std::vector<int>, std::vector<std::vector<acal_match_pair> > > > serializable_nodes;
          for (auto const& item : node_id_map) {
            acal_match_node* node = item.first;
            int node_id = item.second;

            int parent_id;
            if (node->parent_ == 0) {
              // root has a null parent, so give it a special ID
              parent_id = -1;
            } else {
              parent_id = node_id_map[node->parent_];
            }

            std::vector<int> children_ids;
            for (auto const& child_ptr : node->children_) {
              int child_id = node_id_map[child_ptr.get()];
              children_ids.push_back(child_id);
            }
            serializable_nodes[node_id] = std::make_tuple(node->cam_id_, node->node_depth_,
                                                          parent_id, children_ids,
                                                          node->self_to_child_matches_);
          }

          // Use these new node IDs to create serializable representations of each tree
          std::map<size_t, std::map<size_t, std::tuple<size_t, size_t, int> > > serializable_trees;
          for (auto const& item1 : match_trees) {
            std::map<size_t, std::tuple<size_t, size_t, int> > tmp;
            for (auto const& item2 : item1.second) {
              std::shared_ptr<acal_match_tree> tree = item2.second;
              int root_id = node_id_map[tree->root_.get()];
              tmp[item2.first] = std::make_tuple(tree->n_, tree->min_n_tracks_, root_id);
            }
            serializable_trees[item1.first] = tmp;
          }

          // Return a tuple that fully encodes the state of the object
          return py::make_tuple(params, image_paths, all_acams, serializable_vertices,
                                serializable_edges, connected_components_with_ids,
                                focus_tracks, focus_track_metrics, serializable_nodes,
                                serializable_trees, match_tree_metrics);
        },

        [](py::tuple t) {  // __setstate__
          if (t.size() != 11)
            throw std::runtime_error("Can't unpickle acal_match_graph - Invalid state!");

          // Construct vertices
          auto serialized_vertices = t[3].cast<std::map<size_t, std::pair<bool, std::vector<size_t> > > >();
          std::map<size_t, std::shared_ptr<match_vertex> > vertices;
          for (auto const& item : serialized_vertices) {
            size_t vertex_id = item.first;
            auto representation = item.second;

            auto v = std::make_shared<match_vertex>();
            v->cam_id_ = vertex_id;
            v->mark_ = representation.first;
            vertices[vertex_id] = v;
          }

          // Construct edges
          auto serialized_edges = t[4].cast<std::map<size_t, std::tuple<size_t, size_t, std::vector<acal_match_pair> > > >();
          std::map<size_t, std::shared_ptr<match_edge> > edges;
          for (auto const& item : serialized_edges) {
            size_t edge_id = item.first;
            auto representation = item.second;

            auto e = std::make_shared<match_edge>();
            e->id_ = edge_id;
            e->matches_ = std::get<2>(representation);
            edges[edge_id] = e;
          }

          // Construct nodes
          auto serialized_nodes = t[8].cast<std::map<int, std::tuple<size_t, size_t, int, std::vector<int>, std::vector<std::vector<acal_match_pair> > > > >();
          std::map<int, std::shared_ptr<acal_match_node> > nodes;
          for (auto const& item : serialized_nodes) {
            int node_id = item.first;
            auto representation = item.second;

            auto n = std::make_shared<acal_match_node>();
            n->cam_id_ = std::get<0>(representation);
            n->node_depth_ = std::get<1>(representation);
            n->self_to_child_matches_ = std::get<4>(representation);
            nodes[node_id] = n;
          }

          // Construct trees
          auto serialized_trees = t[9].cast<std::map<size_t, std::map<size_t, std::tuple<size_t, size_t, int> > > >();
          std::map<size_t, std::map<size_t, std::shared_ptr<acal_match_tree> > > trees;
          for (auto const& item1 : serialized_trees) {
            std::map<size_t, std::shared_ptr<acal_match_tree> > tmp;
            for (auto const& item2 : item1.second) {
              auto representation = item2.second;

              auto t = std::make_shared<acal_match_tree>();
              t->n_ = std::get<0>(representation);
              t->min_n_tracks_ = std::get<1>(representation);
              tmp[item2.first] = t;
            }
            trees[item1.first] = tmp;
          }

          /* Replace object IDs with object pointers */

          // Put edge pointers into vertices
          for (auto const& item : vertices) {

            size_t vertex_id = item.first;
            const std::shared_ptr<match_vertex>& vertex = item.second;

            std::vector<size_t> edge_ids = serialized_vertices[vertex_id].second;
            for (auto const& edge_id : edge_ids) {
              std::shared_ptr<match_edge>& e = edges[edge_id];
              vertex->add_edge(e.get());
            }
          }

          // Put vertex pointers into edges
          for (auto const& item : edges) {

            size_t edge_id = item.first;
            const std::shared_ptr<match_edge>& edge = item.second;

            size_t v0_id = std::get<0>(serialized_edges[edge_id]);
            size_t v1_id = std::get<1>(serialized_edges[edge_id]);

            const std::shared_ptr<match_vertex>& v0 = vertices[v0_id];
            const std::shared_ptr<match_vertex>& v1 = vertices[v1_id];

            edge->v0_ = v0;
            edge->v1_ = v1;
          }

          // Put node pointers into nodes
          for (auto const& item : nodes) {
            int node_id = item.first;
            const std::shared_ptr<acal_match_node>& node = item.second;

            int parent_id = std::get<2>(serialized_nodes[node_id]);
            std::vector<int> children_ids = std::get<3>(serialized_nodes[node_id]);

            // If this node isn't the root, fill in parent
            if (parent_id != -1) {
              const std::shared_ptr<acal_match_node>& parent = nodes[parent_id];
              node->parent_ = parent.get();
            }

            // Fill in children
            std::vector<std::shared_ptr<acal_match_node> > children;
            for (int child_id : children_ids) {
              const std::shared_ptr<acal_match_node>& child = nodes[child_id];
              children.push_back(child);
            }
            node->children_ = children;
          }

          // Put node pointers into trees
          for (const auto& item1 : trees) {
            for (const auto& item2 : item1.second) {
              int root_id = std::get<2>(serialized_trees[item1.first][item2.first]);
              const std::shared_ptr<acal_match_tree>& tree = item2.second;
              const std::shared_ptr<acal_match_node>& root_node = nodes[root_id];
              tree->root_ = root_node;
            }
          }

          // Put vertex pointers into connected components
          auto connected_components_with_ids = t[5].cast<std::vector<std::vector<size_t> > >();
          std::vector<std::vector<std::shared_ptr<match_vertex> > > connected_components;
          for (auto const& cc_with_ids : connected_components_with_ids) {
            std::vector<std::shared_ptr<match_vertex> > cc;
            for (size_t v_id : cc_with_ids) {
              const std::shared_ptr<match_vertex>& v = vertices[v_id];
              cc.push_back(v);
            }
            connected_components.push_back(cc);
          }

          // Convert the edges map into an edges vector, since that's what acal_match_graph
          // expects. The order of the vector doesn't seem to matter
          std::vector<std::shared_ptr<match_edge> > edges_vector;
          for (auto const& item : edges) {
            const std::shared_ptr<match_edge>& edge = item.second;
            edges_vector.push_back(edge);
          }

          /* Create a new C++ instance */
          acal_match_graph mg;

          /* Set match graph state */
          mg.set_params(t[0].cast<match_params>());
          mg.set_image_paths(t[1].cast<std::map<size_t, std::string> >());
          mg.set_all_acams(t[2].cast<std::map<size_t, vpgl_affine_camera<double> > >());
          mg.set_vertices(vertices);
          mg.set_edges(edges_vector);
          mg.set_connected_components(connected_components);
          mg.set_focus_tracks(t[6].cast<std::map<size_t, std::map<size_t,
                              std::vector< std::map<size_t, vgl_point_2d<double> > > > > >());
          mg.set_focus_track_metrics(t[7].cast<std::vector<double> >());
          mg.set_match_trees(trees);
          mg.set_match_tree_metrics(t[10].cast<std::vector<size_t> >());

          return mg;
        }))
    ;
}

}}  // end pyvxl::acal namespace


PYBIND11_MODULE(_acal, m)
{
  m.doc() =  "Python bindings for the VXL ACAL computer vision library";

  // Need to have access to vgl and vpgl bindings (acal_match_graph uses
  // vpgl_affine_camera and vgl_point_2d)
  py::module importlib_util = py::module::import("importlib.util");
  if (!importlib_util.attr("find_spec")("vxl.vgl").is_none()) {
    py::module::import("vxl.vgl");
  }
  if (!importlib_util.attr("find_spec")("vxl.vpgl").is_none()) {
    py::module::import("vxl.vpgl");
  }

  pyvxl::acal::wrap_f_params(m);
  pyvxl::acal::wrap_match_params(m);
  pyvxl::acal::wrap_acal_corr(m);
  pyvxl::acal::wrap_acal_match_pair(m);
  pyvxl::acal::wrap_match_vertex(m);
  pyvxl::acal::wrap_match_edge(m);
  pyvxl::acal::wrap_acal_match_tree(m);
  pyvxl::acal::wrap_acal_match_graph(m);
}

