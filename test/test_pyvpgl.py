import copy
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

    if hasattr(camA, 'get_matrix'):
      np.testing.assert_allclose(camA.get_matrix(), camB.get_matrix(),
                                rtol = 1e-5, atol = 0)
    else:
      self.assertEqual(camA, camB)


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


# rational camera
class VpglRationalCamera(VpglCameraBase, unittest.TestCase):
  maxDiff = None
  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = vpgl.rational_camera
    self.load = vpgl.load_rational_camera

    # rpc as dictionary
    rpc_dict = {
      "satId": "ABCD",
      "bandId": "P",
      "SpecId": "RPC00B",
      "IMAGE": {
        "errBias": 1.00,
        "errRand": 0.25,
        "lineOffset": 14001,
        "sampOffset": 16001,
        "latOffset": -16.888,
        "longOffset": 179.999,
        "heightOffset": 101,
        "lineScale": 14000,
        "sampScale": 16000,
        "latScale": 0.0783,
        "longScale": 0.0756,
        "heightScale": 100,
        "lineNumCoef": [
          -1.194813E-03, -1.113397E+00, +9.046731E-02, +2.345765E-02,
          +1.821622E-03, -3.800087E-04, -1.804790E-05, +4.634360E-05,
          -4.391746E-04, +6.040494E-06, +5.458677E-06, +1.807921E-06,
          +1.897303E-05, +5.194694E-06, -1.370481E-04, +2.550035E-08,
          -4.761771E-07, -1.230942E-06, -4.227990E-07, -1.129649E-07,
        ],
        "lineDenCoef": [
          +1.000000E+00, -1.438967E-03, +2.176210E-03, -3.936644E-04,
          -9.640095E-05, -1.148574E-06, +5.351234E-06, +1.571059E-05,
          +1.075612E-05, -4.793764E-06, -1.950078E-07, +1.908443E-05,
          +2.007266E-06, +1.937411E-08, -3.958041E-06, -2.356177E-07,
          -1.275287E-08, -9.575094E-07, -1.856964E-08, +0.000000E+00,
        ],
        "sampNumCoef": [
          +8.937241E-03, -3.330377E-03, -1.022227E+00, +1.562245E-02,
          +2.855850E-03, +1.232495E-03, -4.239125E-04, -5.122766E-04,
          -6.330371E-03, -2.184284E-05, +5.048570E-06, -3.106215E-04,
          -1.639877E-04, +1.024656E-06, +2.197754E-04, +1.873784E-05,
          +3.353956E-05, +1.564930E-05, -1.403582E-06, -5.494886E-07,
        ],
        "sampDenCoef": [
          +1.000000E+00, +2.926962E-03, +2.681975E-03, -5.686136E-04,
          -1.195369E-04, -7.231426E-06, +3.135200E-07, +2.028066E-04,
          -5.816487E-06, -3.319875E-05, -6.082773E-07, +1.532235E-06,
          -1.183507E-06, -2.169593E-07, +1.934375E-06, +9.634035E-08,
          +1.260777E-07, +9.890657E-07, +5.879805E-08, +3.260911E-08,
        ],
      },
    }

    # rpc as string
    rpc_string = (
      'satId = "ABCD";\n'
      'bandId = "P";\n'
      'SpecId = "RPC00B";\n'
      'BEGIN_GROUP = IMAGE\n'
      '\terrBias =    1.00;\n'
      '\terrRand =    0.25;\n'
      '\tlineOffset = 14001;\n'
      '\tsampOffset = 16001;\n'
      '\tlatOffset =  -16.888;\n'
      '\tlongOffset =  179.999;\n'
      '\theightOffset = 101;\n'
      '\tlineScale = 14000;\n'
      '\tsampScale = 16000;\n'
      '\tlatScale =    0.0783;\n'
      '\tlongScale =    0.0756;\n'
      '\theightScale = 100;\n'
      '\tlineNumCoef = (\n'
      '\t\t\t-1.194813E-03,\n'
      '\t\t\t-1.113397E+00,\n'
      '\t\t\t+9.046731E-02,\n'
      '\t\t\t+2.345765E-02,\n'
      '\t\t\t+1.821622E-03,\n'
      '\t\t\t-3.800087E-04,\n'
      '\t\t\t-1.804790E-05,\n'
      '\t\t\t+4.634360E-05,\n'
      '\t\t\t-4.391746E-04,\n'
      '\t\t\t+6.040494E-06,\n'
      '\t\t\t+5.458677E-06,\n'
      '\t\t\t+1.807921E-06,\n'
      '\t\t\t+1.897303E-05,\n'
      '\t\t\t+5.194694E-06,\n'
      '\t\t\t-1.370481E-04,\n'
      '\t\t\t+2.550035E-08,\n'
      '\t\t\t-4.761771E-07,\n'
      '\t\t\t-1.230942E-06,\n'
      '\t\t\t-4.227990E-07,\n'
      '\t\t\t-1.129649E-07);\n'
      '\tlineDenCoef = (\n'
      '\t\t\t+1.000000E+00,\n'
      '\t\t\t-1.438967E-03,\n'
      '\t\t\t+2.176210E-03,\n'
      '\t\t\t-3.936644E-04,\n'
      '\t\t\t-9.640095E-05,\n'
      '\t\t\t-1.148574E-06,\n'
      '\t\t\t+5.351234E-06,\n'
      '\t\t\t+1.571059E-05,\n'
      '\t\t\t+1.075612E-05,\n'
      '\t\t\t-4.793764E-06,\n'
      '\t\t\t-1.950078E-07,\n'
      '\t\t\t+1.908443E-05,\n'
      '\t\t\t+2.007266E-06,\n'
      '\t\t\t+1.937411E-08,\n'
      '\t\t\t-3.958041E-06,\n'
      '\t\t\t-2.356177E-07,\n'
      '\t\t\t-1.275287E-08,\n'
      '\t\t\t-9.575094E-07,\n'
      '\t\t\t-1.856964E-08,\n'
      '\t\t\t+0.000000E+00);\n'
      '\tsampNumCoef = (\n'
      '\t\t\t+8.937241E-03,\n'
      '\t\t\t-3.330377E-03,\n'
      '\t\t\t-1.022227E+00,\n'
      '\t\t\t+1.562245E-02,\n'
      '\t\t\t+2.855850E-03,\n'
      '\t\t\t+1.232495E-03,\n'
      '\t\t\t-4.239125E-04,\n'
      '\t\t\t-5.122766E-04,\n'
      '\t\t\t-6.330371E-03,\n'
      '\t\t\t-2.184284E-05,\n'
      '\t\t\t+5.048570E-06,\n'
      '\t\t\t-3.106215E-04,\n'
      '\t\t\t-1.639877E-04,\n'
      '\t\t\t+1.024656E-06,\n'
      '\t\t\t+2.197754E-04,\n'
      '\t\t\t+1.873784E-05,\n'
      '\t\t\t+3.353956E-05,\n'
      '\t\t\t+1.564930E-05,\n'
      '\t\t\t-1.403582E-06,\n'
      '\t\t\t-5.494886E-07);\n'
      '\tsampDenCoef = (\n'
      '\t\t\t+1.000000E+00,\n'
      '\t\t\t+2.926962E-03,\n'
      '\t\t\t+2.681975E-03,\n'
      '\t\t\t-5.686136E-04,\n'
      '\t\t\t-1.195369E-04,\n'
      '\t\t\t-7.231426E-06,\n'
      '\t\t\t+3.135200E-07,\n'
      '\t\t\t+2.028066E-04,\n'
      '\t\t\t-5.816487E-06,\n'
      '\t\t\t-3.319875E-05,\n'
      '\t\t\t-6.082773E-07,\n'
      '\t\t\t+1.532235E-06,\n'
      '\t\t\t-1.183507E-06,\n'
      '\t\t\t-2.169593E-07,\n'
      '\t\t\t+1.934375E-06,\n'
      '\t\t\t+9.634035E-08,\n'
      '\t\t\t+1.260777E-07,\n'
      '\t\t\t+9.890657E-07,\n'
      '\t\t\t+5.879805E-08,\n'
      '\t\t\t+3.260911E-08);\n'
      'END_GROUP = IMAGE\n'
      'END;'
    )

    # test data
    self.test_data = {
      "default": {
        "dict": rpc_dict,
        "string": rpc_string,
      },
    }

  def _create_cam(self, data, check=False):
    cam = vpgl.load_rational_camera_from_str(data['string'])
    if check: self._check_cam(cam, data)
    return cam

  def _check_cam(self, cam, data):

    # result dictionary
    result = copy.deepcopy(data['dict'])
    result.pop('satId')
    result.pop('bandId')
    result['IMAGE'].pop('errBias')
    result['IMAGE'].pop('errRand')

    # convert camera to dictionary
    cam_dict = cam.as_dict(output_rational_order=result['SpecId'])
    cam_dict.pop('satId')
    cam_dict.pop('bandId')

    # compare
    self.assertDictEqual(cam_dict, result)

  def test_project_dict(self):

    # test data
    data = self.test_data['default']

    # RPC origin
    origin_3d = [data['dict']['IMAGE'][key]
                 for key in ('longOffset', 'latOffset', 'heightOffset')]

    # construct cameras from string and dictionary
    camA = vpgl.load_rational_camera_from_str(data['string'])
    self._check_cam(camA, data)
    camB = vpgl.load_rational_camera_from_dict(data['dict'])
    self._check_cam(camB, data)

    # test projection given various offsets from origin
    offsets = [
      (-1e-3, -1e-3, -1),
      (0, 0, 0),
      (1e-3, 1e-3, 1),
    ]
    for offset in offsets:
      pt3d = [origin_3d[i] + offset[i] for i in range(3)]
      ptA = list(camA.project(*pt3d))
      ptB = list(camB.project(*pt3d))
      self.assertEqual(ptA, ptB)


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
