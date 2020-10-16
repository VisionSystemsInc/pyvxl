def load(filename, vil_type):

  # Relative imports, don't pollute vxl.vil import space
  from ._vil import _load_byte, _load_short, _load_float, _load_int
  # from ._vil import _load, _convert_base_sptr_to_byte, _convert_base_sptr_to_short, _convert_base_sptr_to_float, _convert_base_sptr_to_int

  # print("ABOUT TO LOAD FROM PYTHON")
  # file_resource = load_image_resource(filename)
  # base_sptr, native_type = _load(filename)
  # print("DONE LOADING FROM PYTHON")

  # if vil_type is None:
  #   vil_type = native_type

  if vil_type == "byte":
    return _load_byte(filename)
    # return _convert_base_sptr_to_byte(base_sptr)
  elif vil_type == "short":
    return _load_short(filename)
    # return _convert_base_sptr_to_short(base_sptr)
  elif vil_type == "float":
    return _load_float(filename)
    # return _convert_base_sptr_to_float(base_sptr)
  elif vil_type == "int":
    return _load_int(filename)
    # return _convert_base_sptr_to_int(base_sptr)
  else:
    raise ValueError("Unknown vil_type <{}>".format(vil_type))


def stretch_image(image_view, min_limit, max_limit, out_type):

  # Relative imports, don't pollute vxl.vil import space
  from ._vil import _stretch_image_to_byte, _stretch_image_to_short, _stretch_image_to_float

  if (min_limit >= max_limit):
    raise ValueError("vxl.vil.stretch_image: invalid stretch limits")

  # if copy:
  #   # Create new float image
  #   image_view = image_view_float()
  #   image_view.deep_copy(image_view_input)
  # else:
  #   # Modify input view in place, return it
  #   image_view = image_view_input

  if out_type == "byte":
    return _stretch_image_to_byte(image_view, min_limit, max_limit)
  elif out_type == "short":
    return _stretch_image_to_short(image_view, min_limit, max_limit)
  elif out_type == "float":
    return _stretch_image_to_float(image_view, min_limit, max_limit)
  else:
    raise ValueError("vxl.vil.stretch_image: unknown out_type {}".format(out_type))

