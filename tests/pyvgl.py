import unittest

try:
  import numpy as np
except:
  np = None

import vxl

class Vector_2d(unittest.TestCase):

  def test_construct_numbers(self):
    v = vxl.vgl.vector_2d(1.15, -17)

    self.assertEqual(len(v), 2)
    self.assertEqual(v[0], 1.15)
    self.assertEqual(v.x, 1.15)
    self.assertEqual(v[1], -17.0)
    self.assertEqual(v.y, -17.0)

  def test_add(self):
    v = vxl.vgl.vector_2d(1.15, -17)
    u = vxl.vgl.vector_2d(-2.1, 5.2)
    z = v + u

    self.assertEqual(len(z), 2)
    self.assertAlmostEqual(z.x, -0.95)
    self.assertAlmostEqual(z.y, -11.8)

    self.assertEqual(v.x, 1.15)
    self.assertEqual(v.y, -17.0)
    self.assertEqual(u.x, -2.1)
    self.assertEqual(u.y, 5.2)

  def test_subtract(self):
    v = vxl.vgl.vector_2d(1.15, -17)
    u = vxl.vgl.vector_2d(-2.1, 5.2)
    z = v - u

    self.assertEqual(len(z), 2)
    self.assertAlmostEqual(z.x, 3.25)
    self.assertAlmostEqual(z.y, -22.2)

    self.assertEqual(v.x, 1.15)
    self.assertEqual(v.y, -17.0)
    self.assertEqual(u.x, -2.1)
    self.assertEqual(u.y, 5.2)

  @unittest.expectedFailure
  def test_equals(self):
    v = vxl.vgl.vector_2d(1.15, -17)
    u = vxl.vgl.vector_2d(1.15, -17)
    self.assertEqual(v, u)

  @unittest.skipIf(np == None, "Numpy not found")
  def test_construct_numpy(self):
    v = vxl.vgl.vector_2d(np.array([1.15, -17]))

    self.assertEqual(len(v), 2)
    self.assertEqual(v[0], 1.15)
    self.assertEqual(v.x, 1.15)
    self.assertEqual(v[1], -17.0)
    self.assertEqual(v.y, -17.0)

    v = vxl.vgl.vector_2d(np.array([1.15, -17], dtype=int))
    self.assertAlmostEqual(v.x, 1)
    self.assertAlmostEqual(v.y, -17)

class Vector_3d(unittest.TestCase):

  def test_construct_numbers(self):
    v = vxl.vgl.vector_3d(1.15, -17, 42.35)

    self.assertEqual(len(v), 3)
    self.assertEqual(v[0], 1.15)
    self.assertEqual(v.x, 1.15)
    self.assertEqual(v[1], -17.0)
    self.assertEqual(v.y, -17.0)
    self.assertEqual(v[2], 42.35)
    self.assertEqual(v.z, 42.35)

  def test_add(self):
    v = vxl.vgl.vector_3d(1.15, -17, 42.35)
    u = vxl.vgl.vector_3d(-2.1, 5.2, -42.35)
    z = v + u

    self.assertEqual(len(z), 3)
    self.assertAlmostEqual(z.x, -0.95)
    self.assertAlmostEqual(z.y, -11.8)
    self.assertAlmostEqual(z.z, 0)

    self.assertEqual(v.x, 1.15)
    self.assertEqual(v.y, -17.0)
    self.assertEqual(v.z, 42.35)
    self.assertEqual(u.x, -2.1)
    self.assertEqual(u.z, -42.35)

  def test_subtract(self):
    v = vxl.vgl.vector_3d(1.15, -17, 42.35)
    u = vxl.vgl.vector_3d(-2.1, 5.2, -42.35)
    z = v - u

    self.assertEqual(len(z), 3)
    self.assertAlmostEqual(z.x, 3.25)
    self.assertAlmostEqual(z.y, -22.2)
    self.assertAlmostEqual(z.z, 84.7)

    self.assertEqual(v.x, 1.15)
    self.assertEqual(v.y, -17.0)
    self.assertEqual(v.z, 42.35)
    self.assertEqual(u.x, -2.1)
    self.assertEqual(u.z, -42.35)

  @unittest.expectedFailure
  def test_equals(self):
    v = vxl.vgl.vector_3d(1.15, -17, 42.35)
    u = vxl.vgl.vector_3d(-2.1, 5.2, -42.35)
    self.assertEqual(v, u)

  @unittest.skipIf(np == None, "Numpy not found")
  def test_construct_numpy(self):
    v = vxl.vgl.vector_3d(np.array([1.15, -17, 42.35]))

    self.assertEqual(len(v), 3)
    self.assertEqual(v[0], 1.15)
    self.assertEqual(v.x, 1.15)
    self.assertEqual(v[1], -17.0)
    self.assertEqual(v.y, -17.0)
    self.assertEqual(v[2], 42.35)
    self.assertEqual(v.z, 42.35)

    v = vxl.vgl.vector_3d(np.array([1.15, -17, 42.35], dtype=int))
    self.assertAlmostEqual(v.x, 1)
    self.assertAlmostEqual(v.y, -17)
    self.assertAlmostEqual(v.z, 42)

if __name__ == '__main__':
  unittest.main()