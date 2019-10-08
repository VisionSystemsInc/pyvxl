import unittest

try:
  import numpy as np
except:
  np = None

from vxl import vpgl


# class to test expected enumeration values in pybind11 binding
class VpglEnumeration(object):

  def test_enum_values(self):

    # test known enumeration elements
    for item in self.enum_values:
      self.assertIsNotNone(getattr(self.cls, item, None))

    # test invalid enumeration element
    with self.assertRaises(AttributeError):
      self.cls.THIS_IS_AN_INVALID_VALUE


# test vpgl_rational_order enumeration
class VpglRationalOrder(VpglEnumeration, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    self.cls = vpgl.rational_order
    self.enum_values = ['VXL', 'RPC00B']
    super().__init__(*args, **kwargs)


if __name__ == '__main__':
  unittest.main()
