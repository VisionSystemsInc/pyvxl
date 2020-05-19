import numpy as np
import pickle
import unittest

from utils import VxlBase

from vxl import vgl
from vxl import vnl
from vxl import vpgl
from vxl.contrib import acal


class AcalBase(VxlBase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.default_data = dict()

  def test_create(self):
    instance = self.cls()
    self.assertIsInstance(instance, self.cls)
    self.assertAttributes(instance, self.default_data)


class acal_f_params(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.f_params
    self.default_data = {
        'epi_dist_mul': 2.5,
        'max_epi_dist': 5.0,
        'F_similar_abcd_tol': 0.01,
        'F_similar_e_tol': 1.0,
        'ray_uncertainty_tol': 50.0,
        'min_num_matches': 5,
    }


class acal_match_params(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_params
    self.default_data = {
        'min_n_tracks': 3,
        'min_n_cams': 3,
        'max_proj_error': 1.0,
        'max_uncal_proj_error': 20.0,
    }


class acal_corr(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.corr
    self.default_data = {
        'id': np.uint(-1),
        'pt': {'x': -1,
               'y': -1
               }
    }


class acal_match_pair(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_pair
    self.default_data = {
      'corr1': {'id': np.uint(-1),
                'pt': {'x': -1,
                       'y': -1
                       }},
      'corr2': {'id': np.uint(-1),
                'pt': {'x': -1,
                       'y': -1
                       }},
    }


class acal_match_vertex(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_vertex
    self.default_data = {
      'cam_id': np.uint(-1),
      'mark': False,
    }


class acal_match_edge(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_edge
    self.default_data = {
      'id': np.uint(-1),
      'matches': [],
    }


class acal_match_node(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_node
    self.default_data = {
      'cam_id': np.uint(-1),
      'node_depth': 0,
      'children': [],
      'self_to_child_matches': [],
    }


# class acal_match_tree(AcalBase, unittest.TestCase):

#   def __init__(self, *args, **kwargs):
#     super().__init__(*args, **kwargs)
#     self.cls = acal.match_tree
#     import pdb; pdb.set_trace()
#     t = acal.match_tree()
#     del t  # segfaults
#     self.default_data = {
#       'min_n_tracks': 1,
#     }


class acal_match_graph(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_graph


  def construct_example(self):
    incidence_matrix = dict()
    incidence_matrix[0] = dict()
    incidence_matrix[0][1] = list()
    incidence_matrix[0][1].append(acal.match_pair(acal.corr(0, vgl.point_2d(0, 0)), acal.corr(1, vgl.point_2d(1, 1))))
    incidence_matrix[0][1].append(acal.match_pair(acal.corr(0, vgl.point_2d(2, 2)), acal.corr(1, vgl.point_2d(3, 3))))
    incidence_matrix = incidence_matrix

    mg = acal.match_graph(incidence_matrix)
    mg.image_paths = {0: 'image1', 1: 'image2'}

    acam1 = np.full((3, 4), 1)
    acam1[2, :3] = 0
    acam1[0, 0] = 5
    acam2 = np.full((3, 4), 2)
    acam2[2, :3] = 0
    acam2[2, 3] = 1
    acam2[0, 0] = 5
    mg.acams = {0: vpgl.affine_camera(vnl.matrix_fixed_3x4(acam1)),
                1: vpgl.affine_camera(vnl.matrix_fixed_3x4(acam2))}

    mg.find_connected_components()
    mg.compute_focus_tracks()
    mg.compute_match_trees()
    mg.validate_match_trees_and_set_metric()

    return mg


  def test_pickling(self):

    mg = self.construct_example()

    serialized_mg = pickle.dumps(mg)
    reconstructed_mg = pickle.loads(serialized_mg)

    # Make sure vertices were reconstructed correctly
    self.assertTrue(len(reconstructed_mg.vertices) == 2, "Wrong number of vertices")
    self.assertTrue(reconstructed_mg.vertices[0].cam_id == mg.vertices[0].cam_id, "Vertex 0 has the wrong camera ID")
    self.assertTrue(reconstructed_mg.vertices[1].cam_id == mg.vertices[1].cam_id, "Vertex 1 has the wrong camera ID")

    # Make sure edges were reconstructed correctly
    self.assertTrue(len(reconstructed_mg.edges) == 1, "Wrong number of edges")

    self.assertTrue(reconstructed_mg.edges[0].id == mg.edges[0].id, "Edge has the wrong id")

    self.assertTrue(reconstructed_mg.edges[0].v0 == reconstructed_mg.vertices[0], "Edge v0 doesn't match vertex 0")
    self.assertTrue(reconstructed_mg.edges[0].v1 == reconstructed_mg.vertices[1], "Edge v1 doesn't match vertex 1")


    self.assertTrue(reconstructed_mg.edges[0].matches[0].corr1.id == 0, "Edge match 0 corr 1 ID is wrong")
    self.assertTrue(reconstructed_mg.edges[0].matches[0].corr1.pt.x == 0, "Edge match 0 corr 1 x is wrong")
    self.assertTrue(reconstructed_mg.edges[0].matches[0].corr1.pt.y == 0, "Edge match 0 corr 1 y is wrong")

    self.assertTrue(reconstructed_mg.edges[0].matches[0].corr2.id == 1, "Edge match 0 corr 2 ID is wrong")
    self.assertTrue(reconstructed_mg.edges[0].matches[0].corr2.pt.x == 1, "Edge match 0 corr 2 x is wrong")
    self.assertTrue(reconstructed_mg.edges[0].matches[0].corr2.pt.y == 1, "Edge match 0 corr 2 y is wrong")

    self.assertTrue(reconstructed_mg.edges[0].matches[1].corr1.id == 0, "Edge match 1 corr 1 ID is wrong")
    self.assertTrue(reconstructed_mg.edges[0].matches[1].corr1.pt.x == 2, "Edge match 1 corr 1 x is wrong")
    self.assertTrue(reconstructed_mg.edges[0].matches[1].corr1.pt.y == 2, "Edge match 1 corr 1 y is wrong")

    self.assertTrue(reconstructed_mg.edges[0].matches[1].corr2.id == 1, "Edge match 1 corr 2 ID is wrong")
    self.assertTrue(reconstructed_mg.edges[0].matches[1].corr2.pt.x == 3, "Edge match 1 corr 2 x is wrong")
    self.assertTrue(reconstructed_mg.edges[0].matches[1].corr2.pt.y == 3, "Edge match 1 corr 2 y is wrong")

    # Make sure acams were reconstructed correctly
    self.assertTrue(len(reconstructed_mg.acams) == 2, "Wrong number of acams")

    self.assertTrue(np.all(np.array(reconstructed_mg.acams[0].get_matrix()) == np.array([[5, 1, 1, 1],
                                                                                         [1, 1, 1, 1],
                                                                                         [0, 0, 0, 1]])),
                    "Acam 0's matrix is wrong")
    self.assertTrue(np.all(np.array(reconstructed_mg.acams[1].get_matrix()) == np.array([[5, 2, 2, 2],
                                                                                         [2, 2, 2, 2],
                                                                                         [0, 0, 0, 1]])),
                    "Acam 1's matrix is wrong")

    # Make sure connected components were reconstructed correctly
    self.assertTrue(len(reconstructed_mg.connected_components) == 1, "Wrong number of connected components")
    self.assertTrue(len(reconstructed_mg.connected_components[0]) == 2, "Connected component 0 has the wrong size")

    self.assertTrue(reconstructed_mg.connected_components[0][0] == reconstructed_mg.vertices[0],
                    "Connected component 0 has the wrong vertex 0")
    self.assertTrue(reconstructed_mg.connected_components[0][1] == reconstructed_mg.vertices[1],
                    "Connected component 0 has the wrong vertex 1")

    # Make sure trees were reconstructed correctly
    self.assertTrue(len(reconstructed_mg.trees) == 1, "Wrong number of connected components in trees")
    self.assertTrue(len(reconstructed_mg.trees[0]) == 2, "Wrong number of trees in connected component 0")

    self.assertTrue(reconstructed_mg.trees[0][0].root is not None, "Tree 0 is missing its root!")
    self.assertTrue(reconstructed_mg.trees[0][1].root is not None, "Tree 1 is missing its root!")

    self.assertTrue(reconstructed_mg.trees[0][0].root.is_root(), "Tree 0's root isn't a root node")
    self.assertFalse(reconstructed_mg.trees[0][0].root.is_leaf(), "Tree 0's root is a leaf node")

    self.assertTrue(reconstructed_mg.trees[0][1].root.is_root(), "Tree 1's root isn't a root node")
    self.assertFalse(reconstructed_mg.trees[0][1].root.is_leaf(), "Tree 1's root is a leaf node")

    self.assertTrue(reconstructed_mg.trees[0][0].root.cam_id == mg.trees[0][0].root.cam_id,
                    "Tree 0's root has the wrong cam ID")
    self.assertTrue(reconstructed_mg.trees[0][1].root.cam_id == mg.trees[0][1].root.cam_id,
                    "Tree 1's root has the wrong cam ID")

    self.assertTrue(len(reconstructed_mg.trees[0][0].root.children) == 1, "Tree 0 is missing its child!")
    self.assertTrue(len(reconstructed_mg.trees[0][1].root.children) == 1, "Tree 1 is missing its child!")

    self.assertTrue(reconstructed_mg.trees[0][0].root.children[0].cam_id == mg.trees[0][0].root.children[0].cam_id,
                    "Tree 0's child has the wrong cam ID")
    self.assertTrue(reconstructed_mg.trees[0][1].root.children[0].cam_id == mg.trees[0][1].root.children[0].cam_id,
                    "Tree 1's child has the wrong cam ID")

    self.assertFalse(reconstructed_mg.trees[0][0].root.children[0].is_root(),
                     "Tree 0's child is a root node")
    self.assertTrue(reconstructed_mg.trees[0][0].root.children[0].is_leaf(),
                    "Tree 0's child isn't a leaf node")

    self.assertFalse(reconstructed_mg.trees[0][1].root.children[0].is_root(),
                     "Tree 1's child is a root node")
    self.assertTrue(reconstructed_mg.trees[0][1].root.children[0].is_leaf(),
                    "Tree 1's child isn't a leaf node")

    self.assertTrue(reconstructed_mg.trees[0][0].root.children[0].parent(reconstructed_mg.trees[0][0].root) == reconstructed_mg.trees[0][0].root,
                    "Tree 0's child's parent isn't tree 0's root")
    self.assertTrue(reconstructed_mg.trees[0][1].root.children[0].parent(reconstructed_mg.trees[0][1].root) == reconstructed_mg.trees[0][1].root,
                    "Tree 1's child's parent isn't tree 1's root")


if __name__ == '__main__':
  unittest.main()
