import os
import unittest
import utils
import uuid

try:
  import numpy as np
except:
  np = None

from vxl import vil

class VilFile(utils.VxlBase, unittest.TestCase):


  @staticmethod
  def _image_gradient(shape):
    ''' Generate gradient image'''
    bands = list()
    for i in range(shape[2]):
      if i % 2:
        band = np.tile(np.linspace(0, 1, shape[1]), (shape[0], 1))
      else:
        band = np.tile(np.linspace(0, 1, shape[0]), (shape[1], 1)).T
      bands.append(band)
    return np.dstack(bands).squeeze()


  def _test_save_load(self, vil_type='byte', load_type=None, ext='.tif',
                      nbands=1, rng=None, atol=1):
    '''
    vil.save followed by vil.load

    Parameters
    ----------
    self
      This object
    vil_type
      Type of vil object to save
    load_type
      Type of vil object to load, if different from the save type
    ext
      Output extension
    nbands
      Number of output bands/channels
    rng
      Output range ``[min, max]``
    atol
      Absolute tolerance for ``numpy.testing.assert_allclose``
    '''

    # default load type
    if load_type is None:
      load_type = vil_type

    # numpy dtype for save
    dtype_mapping = {
      'byte': np.uint8,
      'short': np.uint16,
      'float': np.float32,
      'int': np.int32,
    }
    dtype = dtype_mapping[vil_type]

    # default range
    if rng is None:
      if np.issubdtype(dtype, np.integer):
        iinfo = np.iinfo(dtype)
        rng = [iinfo.min, iinfo.max]
      else:
        rng = [0.0, 255.0]

    # create floating point image on range [0, 1]
    shape = (256, 256, nbands)
    image = self._image_gradient(shape)

    # stretch to rng
    image = (image * (rng[1] - rng[0])) + rng[0]

    # vil.image_view_*
    vil_func_mapping = {
      'byte': vil.image_view_byte,
      'short': vil.image_view_uint16,
      'float': vil.image_view_float,
      'int': vil.image_view_int,
    }
    vil_func = vil_func_mapping[vil_type]
    view = vil_func(image)

    # save to file
    basename = f'{vil_type}_{nbands}band_{uuid.uuid4()}{ext}'
    file = os.path.join(self.temp_dir.name, basename)

    # save to file
    vil.save_image_view(view, file)
    if not os.path.exists(file):
      raise AssertionError('File was not saved')

    # load from file
    view_load = vil.load(file, load_type)
    image_load = np.array(view_load)

    # check loaded dtype
    _dtype = dtype_mapping[load_type]
    self.assertEqual(image_load.dtype, _dtype,
                     f"vil.load is not expected dtype={_dtype.__name__}")

    # compare loaded values
    np.testing.assert_allclose(image, image_load, rtol=0, atol=atol)


  @unittest.skipUnless(np, 'Numpy not found')
  def test_save_load(self):
    '''Test various vil.save/vil.load combinations'''
    data = [
      # TIFF
      {'ext': '.tif', 'vil_type': 'byte', 'nbands': 1},
      {'ext': '.tif', 'vil_type': 'byte', 'nbands': 3},

      {'ext': '.tif', 'vil_type': 'short', 'nbands': 1},
      {'ext': '.tif', 'vil_type': 'short', 'nbands': 3},

      {'ext': '.tif', 'vil_type': 'float', 'nbands': 1},
      {'ext': '.tif', 'vil_type': 'float', 'nbands': 3},

      # PNG
      {'ext': '.png', 'vil_type': 'byte', 'nbands': 1},
      {'ext': '.png', 'vil_type': 'byte', 'nbands': 3},

      {'ext': '.png', 'vil_type': 'short', 'nbands': 1},
      {'ext': '.png', 'vil_type': 'short', 'nbands': 3},

      # JPEG (increase atol for lossy compression)
      {'ext': '.jpg', 'vil_type': 'byte', 'nbands': 1, 'atol': 5},
      {'ext': '.jpg', 'vil_type': 'byte', 'nbands': 3, 'atol': 5},

      # Save as one type, load as another
      {'vil_type': 'byte', 'load_type': 'short'},
      {'vil_type': 'byte', 'load_type': 'float'},

      {'vil_type': 'short', 'load_type': 'byte', 'rng': [0, 255]},
      {'vil_type': 'short', 'load_type': 'float'},

      {'vil_type': 'float', 'load_type': 'byte', 'rng': [0, 255]},
      {'vil_type': 'float', 'load_type': 'short', 'rng': [0, 65535]},
    ]

    for item in data:
      with self.subTest(item):
        self._test_save_load(**item)


  def test_load_no_file(self):
    '''Raise an exception when loading a file that does not exist'''

    vil_types = ('byte', 'short', 'float', 'int')
    for vil_type in vil_types:
      with self.subTest(vil_type), self.assertRaises(RuntimeError):
        _ = vil.load('/not/a/file', vil_type)
