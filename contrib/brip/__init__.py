from ._brip import *
from vxl import vil


def truncate_nitf_image(img, is_byte, is_scale):
  """Truncate the value of a vil image view to a short or byte.

  Parameters
  ----------
  img : vil_image_view<vxl_uint_16>
    The input image.
  is_byte : bool
    Whether to the returned image will be of type vxl_byte.
    If false, the returned type will be short.
  is_scale : bool
    If is_byte is True, this flag determines whether to scale
    the data values of the returned image to cover the range [0,255].
  """

  from ._brip import _truncate_nitf_image_to_byte, _truncate_nitf_image_to_short

  if is_byte:
    return _truncate_nitf_image_to_byte(img, is_scale)
  else:
    return _truncate_nitf_image_to_short(img)
