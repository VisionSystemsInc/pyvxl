import copy
import numpy as np
import pickle
import unittest
import utils

from vxl import vgl
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


class acal_match_graph(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_graph

  @unittest.skip("not yet implemented")
  def test_equal(self):
    pass


if __name__ == '__main__':
  unittest.main()
