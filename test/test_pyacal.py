import copy
import numpy as np
import pickle
import unittest
import utils

from vxl import vgl
from vxl import vnl
from vxl import vpgl
from vxl.contrib import acal


class AcalBase(utils.VxlBase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.default_data = dict()

  def _set_data_full(self):
    data = copy.deepcopy(self.default_data)
    data = utils.update_nested_dict(data, getattr(self, 'set_data', {}))
    return data

  def test_create(self):
    instance = self.cls()
    self.assertIsInstance(instance, self.cls)
    self.assertAttributes(instance, self.default_data)

  @utils.skipUnlessAttr('set_data')
  def test_init(self):
    init_data = self._set_data_full()
    instance = self.cls(**init_data)
    self.assertIsInstance(instance, self.cls)
    self.assertAttributes(instance, init_data)

  @utils.skipUnlessAttr('set_data')
  def test_set(self):
    instance = self.cls()
    instance.set(**self.set_data)
    new_data = self._set_data_full()
    self.assertAttributes(instance, new_data)

  def test_equal(self):
    init_data = self._set_data_full()
    instance_A = self.cls(**init_data)
    instance_B = self.cls(**init_data)
    self.assertEqual(instance_A, instance_B)

  @utils.skipUnlessClassAttr('__getstate__')
  def test_pickle(self):
    init_data = self._set_data_full()
    insstance_A = self.cls(**init_data)
    insstance_B = pickle.loads(pickle.dumps(insstance_A))
    self.assertEqual(insstance_A, insstance_B)


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
    self.set_data = {
        'epi_dist_mul': 3.0,
        'max_epi_dist': 5.5,
        'F_similar_abcd_tol': 0.51,
        'F_similar_e_tol': 1.5,
        'ray_uncertainty_tol': 50.5,
        'min_num_matches': 6,
    }


class acal_corr(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.corr
    self.default_data = {
        'id': np.uint(-1),
        'pt': vgl.point_2d(-1, -1)
    }
    self.set_data = {
        'id': 10,
        'pt': vgl.point_2d(10.0, 20.0)
    }


class acal_match_pair(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_pair
    self.default_data = {
        'corr1': {'id': np.uint(-1), 'pt': vgl.point_2d(-1, -1)},
        'corr2': {'id': np.uint(-1), 'pt': vgl.point_2d(-1, -1)},
    }
    self.set_data = {
        'corr1': {'id': 10, 'pt': vgl.point_2d(10.0, 20.0)},
        'corr2': {'id': 15, 'pt': vgl.point_2d(15.0, 25.0)},
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
    self.set_data = {
        'min_n_tracks': 5,
        'min_n_cams': 10,
        'max_proj_error': 10.0,
        'max_uncal_proj_error': 30.0,
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


class acal_match_tree(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_tree
    self.default_data = {
      'min_n_tracks': 1,
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


class acal_match_graph(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_graph

  @unittest.skip("not yet implemented")
  def test_equal(self):
    pass


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


if __name__ == '__main__':
  unittest.main()
