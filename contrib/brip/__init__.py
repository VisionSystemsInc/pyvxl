from ._brip import *

def truncate_nitf_image(img, is_byte, is_scale):
  if is_byte:
    return truncate_nitf_image_to_byte(img, is_scale)
  else:
    return truncate_nitf_image_to_short(img)
