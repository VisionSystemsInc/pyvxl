import unittest
import os
import pickle
import tempfile

try:
  import numpy as np
except:
  np = None

from vxl import vgl
import vxl.vgl.algo
from vxl import vnl
from vxl import vpgl


# helper function: populate vnl.matrix_fixed
# (without numpy requirement - input is list of lists)
def _matrix_fixed(data, vnl_cls):

  matrix = vnl_cls()
  sz = matrix.shape

  if len(data) != sz[0] or any(len(row) != sz[1] for row in data):
    raise Exception('fixed matrix expects {}x{} data'.format(*sz))

  for r in range(sz[0]):
    for c in range(sz[1]):
      matrix[r, c] = data[r][c]

  return matrix

def matrix_fixed_3x3(data):
  return _matrix_fixed(data, vnl.matrix_fixed_3x3)

def matrix_fixed_3x4(data):
  return _matrix_fixed(data, vnl.matrix_fixed_3x4)


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

  @unittest.skipUnless(np, "Numpy not found")
  def test_save_load(self):
    # save/load to an ASCII file of limited precision
    # test near equality of camera matrices via "assert_allclose"

    data = self.test_data['default']
    camA = self._create_cam(data, check = True)

    with tempfile.TemporaryDirectory() as temp_dir:
      file = os.path.join(temp_dir, 'camera.txt')
      camA.save(file)
      camB = self.load(file)

    np.testing.assert_allclose(camA.get_matrix(), camB.get_matrix(),
                               rtol = 1e-5, atol = 0)


# projective camera
class VpglProjCamera(VpglCameraBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = vpgl.proj_camera
    self.load = vpgl.load_proj_camera

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
    self.load = vpgl.load_affine_camera

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


# perspective camera
class VpglPerspectiveCamera(VpglCameraBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = vpgl.perspective_camera
    self.load = vpgl.load_perspective_camera

    # test data
    self.test_data = {
      "default": {
        "calibration": vpgl.calibration_matrix(matrix_fixed_3x3([
            [2000.0, 0.0, 512.0],
            [0.0, 2000.0, 384.0],
            [0.0, 0.0, 1.0]
        ])),
        "camera_center": vgl.point_3d(0.0, 0.0, -10.0),
        "rotation": vgl.algo.rotation_3d(),
      }
    }

  def _create_cam(self, data, check = False):
    cam = self.cls(data['calibration'], data['camera_center'], data['rotation'])
    if check: self._check_cam(cam, data)
    return cam

  def _check_cam(self, cam, data):
    self.assertEqual(cam.calibration, data['calibration'])
    self.assertEqual(cam.camera_center, data['camera_center'])
    self.assertEqual(cam.rotation, data['rotation'])


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
