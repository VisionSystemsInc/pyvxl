import unittest
import pickle

try:
  import numpy as np
except:
  np = None

from vxl import vnl


# ----------
# VNL VECTOR
# ----------
class VnlVectorBase(object):

  def __init__(self, *args, **kwargs):
    self.fixed_size = None
    super().__init__(*args, **kwargs)

  @unittest.skipUnless(np, "Numpy not found")
  def test_construct_numpy(self):
    data = [2.51, 3.14, -7.8, 0, 10.0]
    data = np.array(data, dtype=np.double)

    size = self.fixed_size or data.size
    data = data[0:size]

    x = self.cls(data)

    self.assertEqual(len(x), size)
    self.assertEqual(x.size, size)

    for i in range(size):
      self.assertEqual(x[i], data[i])
      self.assertEqual(x.get(i), data[i])

    for i in range(-size, 0):
      self.assertEqual(x[i], data[i])

  def test_construct_len(self):
    size = self.fixed_size or 5

    if self.fixed_size is None:
      x = self.cls(size)
    else:
      x = self.cls()

    self.assertEqual(len(x), size)
    self.assertEqual(x.size, size)

  def test_construct_len_value(self):
    size = self.fixed_size or 6
    value = 1.23

    if self.fixed_size is None:
      x = self.cls(size, value)
    else:
      x = self.cls(value)

    self.assertEqual(len(x), size)
    self.assertEqual(x.size, size)

    for i in range(size):
      self.assertEqual(x[i], value)
      self.assertEqual(x.get(i),value)

  def test_setitem(self):
    size = self.fixed_size or 5

    if self.fixed_size is None:
      x = self.cls(size)
    else:
      x = self.cls()

    for i in range(size):
      x[i] = i
      self.assertEqual(x[i], i)

  def test_operator_equal(self):
    size = self.fixed_size or 4
    value = 1.23

    if self.fixed_size is None:
      a = self.cls(size, value)
      b = self.cls(size, value)
    else:
      a = self.cls(value)
      b = self.cls(value)

    self.assertEqual(a,b)

  def test_pickle(self):
    size = self.fixed_size or 4
    value = 1.23

    if self.fixed_size is None:
      a = self.cls(size, value)
    else:
      a = self.cls(value)

    b = pickle.loads(pickle.dumps(a))
    self.assertEqual(a,b)


class VnlVector(VnlVectorBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vnl.vector
    super().__init__(*args, **kwargs)

class VnlVectorFixed3(VnlVectorBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vnl.vector_fixed_3
    super().__init__(*args, **kwargs)
    self.fixed_size = 3

class VnlVectorFixed4(VnlVectorBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vnl.vector_fixed_4
    super().__init__(*args, **kwargs)
    self.fixed_size = 4


# ----------
# VNL MATRIX
# ----------
class VnlMatrixBase(object):

  def __init__(self, *args, **kwargs):
    self.fixed_shape = None
    super().__init__(*args, **kwargs)

  @unittest.skipUnless(np, "Numpy not found")
  def test_construct_numpy(self):
    data = np.zeros((20,20), dtype = np.double)
    data[0:5, 0:5] = [[0.0, 0.1, 0.2, 0.3, 0.4],
                      [1.0, 1.1, 1.2, 1.3, 1.4],
                      [2.0, 2.1, 2.2, 2.3, 2.4],
                      [3.0, 3.1, 3.2, 3.3, 3.4],
                      [4.0, 4.1, 4.2, 4.3, 4.4]]

    shape = self.fixed_shape or data.shape
    data = data[0:shape[0], 0:shape[1]]
    x = self.cls(data)

    self.assertEqual(len(x), shape[0])
    self.assertEqual(x.shape, shape)

    for r in range(shape[0]):
      for c in range(shape[1]):
        self.assertEqual(x[r,c], data[r,c])
        self.assertEqual(x[r][c], data[r,c])
        self.assertEqual(x.get(r,c), data[r,c])

    for r in range(-shape[0], 0):
      for c in range(-shape[1], 0):
        self.assertEqual(x[r,c], data[r,c])

  def test_construct_shape(self):
    shape = self.fixed_shape or (4,5)

    if self.fixed_shape is None:
      x = self.cls(*shape)
    else:
      x = self.cls()

    self.assertEqual(len(x), shape[0])
    self.assertEqual(x.shape, shape)

  def test_construct_shape_value(self):
    shape = self.fixed_shape or (5,6)
    value = 2.345

    if self.fixed_shape is None:
      x = self.cls(*shape, value)
    else:
      x = self.cls(value)

    self.assertEqual(len(x), shape[0])
    self.assertEqual(x.shape, shape)

    for r in range(shape[0]):
      for c in range(shape[1]):
        self.assertEqual(x[r,c], value)
        self.assertEqual(x[r][c], value)
        self.assertEqual(x.get(r,c), value)

  def test_setitem(self):
    shape = self.fixed_shape or (4,5)

    if self.fixed_shape is None:
      x = self.cls(*shape)
    else:
      x = self.cls()

    for r in range(shape[0]):
      for c in range(shape[1]):
        value = shape[0]*r + c
        x[r,c] = value
        self.assertEqual(x[r,c], value)

  def test_operator_equal(self):
    shape = self.fixed_shape or (6,7)
    value = 2.345

    if self.fixed_shape is None:
      a = self.cls(*shape, value)
      b = self.cls(*shape, value)
    else:
      a = self.cls(value)
      b = self.cls(value)

    self.assertEqual(a,b)

  def test_pickle(self):
    shape = self.fixed_shape or (6,7)
    value = 2.345

    if self.fixed_shape is None:
      a = self.cls(*shape, value)
    else:
      a = self.cls(value)

    b = pickle.loads(pickle.dumps(a))
    self.assertEqual(a,b)


class VnlMatrix(VnlMatrixBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vnl.matrix
    super().__init__(*args, **kwargs)

class VnlMatrix3x3(VnlMatrixBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vnl.matrix_fixed_3x3
    super().__init__(*args, **kwargs)
    self.fixed_shape = (3,3)

class VnlMatrix3x4(VnlMatrixBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vnl.matrix_fixed_3x4
    super().__init__(*args, **kwargs)
    self.fixed_shape = (3,4)

class VnlMatrix4x20(VnlMatrixBase, unittest.TestCase):
  def __init__(self, *args, **kwargs):
    self.cls = vnl.matrix_fixed_4x20
    super().__init__(*args, **kwargs)
    self.fixed_shape = (4,20)

if __name__ == '__main__':
  unittest.main()
