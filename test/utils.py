import math

# safe nan check
def safe_isnan(val):
  try:
    return math.isnan(val)
  except TypeError:
    return False

# generic base unittest object
class VxlBase(object):

  # check for attribute
  def assertHasAttr(self, obj, attr, msg = None):
    self.assertTrue(hasattr(obj, attr), msg)

  # check for NaN
  def assertNan(self, expr, msg = None):
    self.assertTrue(safe_isnan(expr), msg)

  # recursive check for attribute existence/value
  def assertAttributes(self, instance, attr_dct):

    def nested_func(obj, dct, parent_keys = tuple()):

      for key, val in dct.items():
        keys = parent_keys + (key,)
        keys_str = '.'.join(keys)

        # check for key existence
        self.assertHasAttr(obj, key,
                           msg = 'attribute {} not found'.format(keys_str))

        # get key value
        obj_val = getattr(obj, key)

        # recursion
        if isinstance(val, dict):
          nested_func(obj_val, val, keys)
          continue

        # check value
        if safe_isnan(val):
          assertFunc = self.assertNan
          opts = {}

        elif isinstance(val, float):
          assertFunc = self.assertAlmostEqual
          opts = {'second': val, 'places': 7}

        else:
          assertFunc = self.assertEqual
          opts = {'second': val}

        opts['msg'] = 'attribute {} unexpected value'.format(keys_str)
        assertFunc(obj_val, **opts)

    # call recursive check
    nested_func(instance, attr_dct)

