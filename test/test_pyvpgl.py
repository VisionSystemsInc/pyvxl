import unittest
import pickle

from vxl import vnl
from vxl import vpgl

# helper function: populate vnl.matrix_fixed_3x4
# (without numpy requirement - input is list of lists)
def matrix_fixed_3x4(data):

  # check size
  if len(data) != 3 or any(len(row) != 4 for row in data):
    raise Exception('3x4 matrix expects 3x4 data')

  matrix = vnl.matrix_fixed_3x4()
  for r in range(matrix.shape[0]):
    for c in range(matrix.shape[1]):
      matrix[r,c] = data[r][c]

  return matrix


# ----------
# VPGL CAMERA
# ----------

# generic camera test class
class VpglCameraBase(object):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)

  def test_create(self):
    data = self.test_data['default']
    cam = self._create_cam(data)
    self._check_cam(cam, data)

  def test_equal(self):
    data = self.test_data['default']
    camA = self._create_cam(data, check = True)
    camB = self._create_cam(data, check = True)
    self.assertEqual(camA, camB)

  def test_clone(self):
    data = self.test_data['default']
    camA = self._create_cam(data, check = True)
    camB = camA
    self.assertIs(camA, camB)
    camC = camA.clone()
    self.assertIsNot(camA, camC)

  def test_pickle(self):
    data = self.test_data['default']
    camA = self._create_cam(data, check = True)
    camB = pickle.loads(pickle.dumps(camA))
    self.assertEqual(camA, camB)

# projective camera
class VpglProjCamera(VpglCameraBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = vpgl.proj_camera

    # test data
    self.test_data = {
      "default": {
        "matrix": matrix_fixed_3x4([
            [437.5, 128.5575, -153.20889, 20153.20898],
            [0.0,  -206.5869, -434.42847, 20434.42968],
            [0.0,   0.642787,   -0.76604,    100.7660]
          ]),
        }
      }

  def _create_cam(self, data, check = False):
    cam = self.cls(data['matrix'])
    if check: self._check_cam(cam, data)
    return cam

  def _check_cam(self, cam, data):
    self.assertEqual(cam.get_matrix(), data['matrix'])

# affine camera
class VpglAffineCamera(VpglCameraBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = vpgl.affine_camera

    # test data
    self.test_data = {
      "default": {
        "matrix": matrix_fixed_3x4([
            [1.34832357134959380,  0.0038174980872772743,  0.27647870881886161,   8.8599950663932052],
            [0.21806927892797245, -0.9263109114580021500, -1.00105353309762050, 538.9337651820003400],
            [0.0, 0.0, 0.0, 1.0]
          ]),
        "viewing_distance": 9325.6025071654913,
        },
      }

  def _create_cam(self, data, check = False):
    cam = self.cls(data['matrix'])
    cam.viewing_distance = data['viewing_distance']
    if check: self._check_cam(cam, data)
    return cam

  def _check_cam(self, cam, data):
    self.assertEqual(cam.get_matrix(), data['matrix'])
    self.assertEqual(cam.viewing_distance, data['viewing_distance'])


# ----------
# OTHER TESTS
# ----------

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
