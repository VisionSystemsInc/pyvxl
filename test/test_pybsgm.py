import unittest

from vxl import vgl
from vxl.contrib import bsgm


class bsgm_generic(object):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.default_data = dict()


  # recursive check for attribute existence/value
  def assertAttributes(self, instance, attr_dct):

    def nested_func(obj, dct, parent_keys = tuple()):

      for key, val in dct.items():
        keys = parent_keys + (key,)
        keys_str = '.'.join(keys)

        # check for key existence
        self.assertTrue(hasattr(obj, key),
                        msg = 'attribute {} not found'.format(keys_str))

        # get key value
        obj_val = getattr(obj, key)

        # recursion
        if isinstance(val, dict):
          nested_func(obj_val, val, keys)
          continue

        # check value
        if isinstance(val, float):
          assertFunc = self.assertAlmostEqual
          opts = {'places': 7}
        else:
          assertFunc = self.assertEqual
          opts = {}

        assertFunc(obj_val, val,
                   msg = 'attribute {} unexpected value'.format(keys_str),
                   **opts)

    # call recursive check
    nested_func(instance, attr_dct)


  def test_create(self):
    instance = self.cls()
    self.assertIsInstance(instance, self.cls)
    self.assertAttributes(instance, self.default_data)


class bsgm_disparity_estimator_params(bsgm_generic, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = bsgm.disparity_estimator_params

    self.default_data = {
        'use_16_directions': False,
        'p1_scale': 1.0,
        'p2_scale': 1.0,
        'use_gradient_weighted_smoothing': True,
        'max_grad': 32.0,
        'perform_quadratic_interp': True,
        'error_check_mode': 1,
        'shadow_thresh': 0,
        'bias_weight': 0.0,
        'bias_dir': vgl.vector_2d_float(1.0, 0.0),
        'census_weight': 0.3,
        'xgrad_weight': 0.7,
        'census_tol': 2,
        'census_rad': 2,
        'print_timing': False,
    }


class bsgm_pairwise_params(bsgm_generic, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = bsgm.pairwise_params

    self.default_data = {
        'shadow_thresh': 20,
        'quad_interp': False,
        'multi_scale_mode': 0,
        'active_disparity_factor': 0.5,
        'downscale_exponent': 2,
        'ground_sample_dist': 0.3,
        'upsample_scale_factor': 1.0,
        'std_dev': 3.75 * 0.3,
        'use_z_vs_d_prob': False,
        'min_z_vs_d_scale': 1.0,
        'z_vs_d_std_dev': 1.0,
        'min_neighbors': 3,
        'max_neighbors': 5,
        'neighbor_dist_factor': 3.0,
    }

    self.default_data['de_params'] = {
        'use_16_directions': False,
        'p1_scale': 1.0,
        'p2_scale': 1.0,
        'use_gradient_weighted_smoothing': True,
        'max_grad': 32.0,
        'perform_quadratic_interp': self.default_data['quad_interp'],
        'error_check_mode': 1,
        'shadow_thresh': self.default_data['shadow_thresh'],
        'bias_weight': 0.0,
        'bias_dir': vgl.vector_2d_float(1.0, 0.0),
        'census_weight': 0.3,
        'xgrad_weight': 0.7,
        'census_tol': 2,
        'census_rad': 2,
        'print_timing': False,
    }


if __name__ == '__main__':
  unittest.main()
