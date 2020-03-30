import collections.abc
import math

# nested dictionary update
def update_nested_dict(d, u):
  for k, v in u.items():
    if isinstance(v, collections.abc.Mapping):
      d[k] = update_nested_dict(d.get(k, {}), v)
    else:
      d[k] = v
  return d

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

        # recursion ("val" dictionary cannot be empty)
        if isinstance(val, dict):
          if not val:
            raise ValueError('empty attribute validator {}'.format(keys_str))
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


# decorator: skip unit test if attribute is missing or None
def skipUnlessAttr(attr):
  def decorator(func):
    def wrapper(self, *args, **kwargs):
      if getattr(self, attr, None) is None:
        self.skipTest('"{}" not available'.format(attr))
      else:
        func(self, *args, **kwargs)
    return wrapper
  return decorator
