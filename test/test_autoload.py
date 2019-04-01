import unittest

try:
  import numpy as np
except:
  np = None


class TestAutoLoad(unittest.TestCase):

  def test_imports(self):
    import vxl
    vxl.vgl
    vxl.vgl.algo
    vxl.vil
    vxl.vnl
    vxl.vpgl
    vxl.vpgl.algo
    vxl.contrib


if __name__ == '__main__':
  unittest.main()
