import unittest

try:
  import numpy as np
except:
  np = None

import vxl

class Vpgl:
  def test_ray_casting(self):

      P = np.array([7.0931,  -2.9319,   0.2307, 173.85],
          [0.7152,  -0.0537, -10.2546, 160.956],
          [0.0041,   0.0045,  -0.0006,   0.3261])

      pcam = vxl.vpgl.proj_camera(P)

      l3d = pcam.backproject(p)

      pass

if __name__ == '__main__':
  unittest.main()