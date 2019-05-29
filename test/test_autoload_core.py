import unittest


class TestAutoLoad(unittest.TestCase):
  """
  Test that importing vxl auto-loads all of the
  core vxl modules.
  """

  def test_imports(self):
    import vxl
    vxl.vgl
    vxl.vgl.algo
    vxl.vil
    vxl.vnl
    vxl.vpgl
    vxl.vpgl.algo


if __name__ == '__main__':
  unittest.main()
