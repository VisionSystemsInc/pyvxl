import unittest

try:
  import numpy as np
except:
  np = None

import vxl

class VglBase(object):

  def test_construct_numbers(self):
    if self.length == 2:
      v = self.cls(1.15, -17)
    else:
      v = self.cls(1.15, -17, 42.35)

    self.assertEqual(len(v), self.length)
    self.assertEqual(v[0], 1.15)
    self.assertEqual(v.x, 1.15)
    self.assertEqual(v[1], -17.0)
    self.assertEqual(v.y, -17.0)
    if self.length == 3:
      self.assertEqual(v[2], 42.35)
      self.assertEqual(v.z, 42.35)

  def test_add(self):
    if self.length == 2:
      v = self.cls(1.15, -17)
      u = self.cls(-2.1, 5.2)
    else:
      v = self.cls(1.15, -17, 42.35)
      u = self.cls(-2.1, 5.2, -42.35)
    z = v + u

    self.assertEqual(len(z), self.length)
    self.assertAlmostEqual(z.x, -0.95)
    self.assertAlmostEqual(z.y, -11.8)
    if self.length == 3:
      self.assertAlmostEqual(z.z, 0)

    self.assertEqual(v.x, 1.15)
    self.assertEqual(u.x, -2.1)
    self.assertEqual(v.y, -17.0)
    self.assertEqual(u.y, 5.2)
    if self.length == 3:
      self.assertEqual(v.z, 42.35)
      self.assertEqual(u.z, -42.35)

  def test_subtract(self):
    if self.length == 2:
      v = self.cls(1.15, -17)
      u = self.cls(-2.1, 5.2)
    else:
      v = self.cls(1.15, -17, 42.35)
      u = self.cls(-2.1, 5.2, -42.35)
    z = v - u

    self.assertEqual(len(z), self.length)
    self.assertAlmostEqual(z.x, 3.25)
    self.assertAlmostEqual(z.y, -22.2)
    if self.length == 3:
      self.assertAlmostEqual(z.z, 84.7)

    self.assertEqual(v.x, 1.15)
    self.assertEqual(u.x, -2.1)
    self.assertEqual(v.y, -17.0)
    self.assertEqual(u.y, 5.2)
    if self.length == 3:
      self.assertEqual(v.z, 42.35)
      self.assertEqual(u.z, -42.35)

  @unittest.expectedFailure
  def test_equals(self):
    if self.length == 2:
      v = self.cls(1.15, -17)
      u = self.cls(1.15, -17)
    else:
      v = self.cls(1.15, -17, 42.35)
      u = self.cls(1.15, -17, 42.35)
    self.assertEqual(v, u)

  @unittest.skipUnless(np, "Numpy not found")
  def test_construct_numpy(self):
    if self.length == 2:
      v = self.cls(np.array([1.15, -17]))
      u = self.cls(np.array([1.15, -17], dtype=int))
    else:
      v = self.cls(np.array([1.15, -17, 42.35]))
      u = self.cls(np.array([1.15, -17, 42.35], dtype=int))

    self.assertEqual(len(v), self.length)
    self.assertEqual(v[0], 1.15)
    self.assertEqual(v.x, 1.15)
    self.assertEqual(v[1], -17.0)
    self.assertEqual(v.y, -17.0)
    if self.length == 3:
      self.assertEqual(v[2], 42.35)
      self.assertEqual(v.z, 42.35)

    self.assertAlmostEqual(u.x, 1)
    self.assertAlmostEqual(u.y, -17)
    if self.length == 3:
      self.assertAlmostEqual(u.z, 42)

class Vector_2d(VglBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vxl.vgl.vector_2d
    self.length = 2
    super().__init__(*args, **kwargs)

class Vector_3d(VglBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vxl.vgl.vector_3d
    self.length = 3
    super().__init__(*args, **kwargs)

class Point_2d(VglBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vxl.vgl.point_2d
    self.length = 2
    super().__init__(*args, **kwargs)

  @unittest.skip("Not implemented")
  def test_add(self):
    pass

class Point_3d(VglBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vxl.vgl.point_3d
    self.length = 3
    super().__init__(*args, **kwargs)

  @unittest.skip("Not implemented")
  def test_add(self):
    pass

if __name__ == '__main__':
  unittest.main()