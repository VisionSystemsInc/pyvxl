import unittest


class TestAutoLoad(unittest.TestCase):
  """
  Test that importing vxl auto-loads all of the
  contrib modules (if they were built).
  """

  def test_imports(self):
    import vxl
    vxl.contrib
    vxl.contrib.bpgl
    vxl.contrib.bpgl.algo
    vxl.contrib.brad
    vxl.contrib.brip
    vxl.contrib.bvxm
    vxl.contrib.bvxm.algo
    vxl.contrib.sdet
    vxl.contrib.sdet.algo


if __name__ == '__main__':
  unittest.main()
