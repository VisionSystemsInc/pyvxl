from ._brip import *


def truncate_nitf_image(img, is_byte, is_scale):
  """exp

  Parameters
  ----------
  img : vil_image_view<vxl_uint_16>
    The input image.
  is_byte : bool
    Whether to the returned image will be of type vxl_byte
  is_scale : bool
    Whether to scale the data values of the returned image to lie between 0 and 255.
  """
  if is_byte:
    return truncate_nitf_image_to_byte(img, is_scale)
  else:
    return truncate_nitf_image_to_short(img)
