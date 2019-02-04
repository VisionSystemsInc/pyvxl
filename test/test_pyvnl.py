import unittest

try:
  import numpy as np
except:
  np = None

import vxl

class VnlVectorBase(object):
  @unittest.skipUnless(np, "Numpy not found")
  def test_construct_numpy(self):
    x = self.cls(np.array([2.51, 3.14, -7.8, 0], dtype=np.double))

    self.assertEqual(len(x), 4)
    self.assertEqual(x.size(), 4)
    self.assertEqual(x[0], 2.51)
    self.assertEqual(x.get(0), 2.51)
    self.assertEqual(x[1], 3.14)
    self.assertEqual(x.get(1), 3.14)
    self.assertEqual(x[2], -7.8)
    self.assertEqual(x.get(2), -7.8)
    self.assertEqual(x[3], 0.0)
    self.assertEqual(x.get(3), 0.0)

    self.assertEqual(x[-4], 2.51)
    self.assertEqual(x[-3], 3.14)
    self.assertEqual(x[-2], -7.8)
    self.assertEqual(x[-1], 0.0)

  def test_construct_len(self):
    y = self.cls(5)

    self.assertEqual(len(y), 5)
    self.assertEqual(y.size(), 5)

  def test_construct_len_value(self):
    z = self.cls(6,1.23)

    self.assertEqual(len(z), 6)
    self.assertEqual(z.size(), 6)
    for x in range(6):
      self.assertEqual(z[x], 1.23)
      self.assertEqual(z.get(x), 1.23)

class VnlVector(VnlVectorBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vxl.vnl.vector
    super().__init__(*args, **kwargs)

if __name__ == '__main__':
  unittest.main()