import unittest

from vxl.contrib import acal


class acal_generic(object):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.default_data = dict()

  def assertHasAttr(self, obj, attr):
    self.assertTrue(hasattr(obj, attr),
                    msg = 'attribute {} not found'.format(attr))

  def test_create(self):
    instance = self.cls()
    self.assertIsInstance(instance, self.cls)

    for k, v in self.default_data.items():
      self.assertHasAttr(instance, k)
      self.assertEqual(getattr(instance, k), v,
                       '{} has unexpected value'.format(k))


class acal_f_params(acal_generic, unittest.TestCase):

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


class acal_match_params(acal_generic, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_params
    self.default_data = {
        'min_n_tracks': 3,
        'min_n_cams': 3,
        'max_proj_error': 1.0,
        'max_uncal_proj_error': 20.0,
    }


class acal_match_graph(acal_generic, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_graph


if __name__ == '__main__':
  unittest.main()
