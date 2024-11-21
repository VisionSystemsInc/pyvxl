from ._vpgl import *
from vxl import vgl, vnl


def load_rational_camera(cam_fname):
  """
  Loads either a local rational camera, or a rational camera, whichever is appropriate.

  Parameters
  ----------
  cam_fname : string
    The file name of the camera file (e.g. .rpb).

  Returns
  -------
  vpgl_rational_camera

  Raises
  ------
  ValueError
    If cam_fname is not a valid camera file.
  """
  from ._vpgl import _read_local_rational_camera, _read_rational_camera

  local_cam = _read_local_rational_camera(cam_fname)
  if (local_cam is not None):
    return local_cam
  else:
    cam = _read_rational_camera(cam_fname)
    if (cam is not None):
      return cam
    else:
      raise ValueError("Invalid camera file <{}>".format(cam_fname))


def load_rational_camera_from_txt(cam_fname):
  """
  Loads a rational camera.

  Parameters
  ----------
  cam_fname : string
    The file name of the camera file.

  Returns
  -------
  vpgl_rational_camera

  Raises
  ------
  ValueError
    If cam_fname is not a valid camera file.
  """

  from ._vpgl import _read_rational_camera_from_txt

  cam = _read_rational_camera_from_txt(cam_fname)
  if cam:
    return cam
  else:
    raise ValueError("Invalid camera file <{}>".format(cam_fname))


def load_rational_camera_from_dict(rpc_dict):
  """
  Load rational camera from a python dictionary

  Parameters
  ----------
  rpc_dict : :obj:`dict`
    Dictionary of RPC parameters
  """

  # coefficient order
  spec_id = rpc_dict["SpecId"]
  input_rational_order = getattr(_vpgl.rational_order, spec_id.upper(), None)
  if not input_rational_order:
    raise RuntimeError(f"Invalid rational camera ordering {spec_id}")

  # create VXL object
  IMAGE = rpc_dict["IMAGE"]
  camera = _vpgl.rational_camera(

    # rational camera ordering
    input_rational_order=input_rational_order,

    # offset
    u_off=IMAGE["sampOffset"],
    v_off=IMAGE["lineOffset"],
    x_off=IMAGE["longOffset"],
    y_off=IMAGE["latOffset"],
    z_off=IMAGE["heightOffset"],

    # scale
    u_scale=IMAGE["sampScale"],
    v_scale=IMAGE["lineScale"],
    x_scale=IMAGE["longScale"],
    y_scale=IMAGE["latScale"],
    z_scale=IMAGE["heightScale"],

    # coefficients
    neu_u=IMAGE["sampNumCoef"],
    den_u=IMAGE["sampDenCoef"],
    neu_v=IMAGE["lineNumCoef"],
    den_v=IMAGE["lineDenCoef"],

  )

  return camera


def rational_camera_as_dict(cam, output_rational_order='RPC00B'):
  """
  Convert rational camera to python dictionary
  """

  # coefficient order
  if isinstance(output_rational_order, str):
    output_rational_order = getattr(_vpgl.rational_order, output_rational_order.upper())

  # coefficients
  coefficient_matrix = cam.coefficient_matrix(output_rational_order)
  poly_mapping = {
    "sampNumCoef": _vpgl.rational_camera.poly_index.NEU_U,
    "sampDenCoef": _vpgl.rational_camera.poly_index.DEN_U,
    "lineNumCoef": _vpgl.rational_camera.poly_index.NEU_V,
    "lineDenCoef": _vpgl.rational_camera.poly_index.DEN_V,
  }

  def _coefficients(name):
    return list(coefficient_matrix[poly_mapping[name]])

  # offset & scale mapping
  coord_mapping = {
    "samp":   _vpgl.rational_camera.coor_index.U_INDX,
    "line":   _vpgl.rational_camera.coor_index.V_INDX,
    "long":   _vpgl.rational_camera.coor_index.X_INDX,
    "lat":    _vpgl.rational_camera.coor_index.Y_INDX,
    "height": _vpgl.rational_camera.coor_index.Z_INDX,
  }

  def _offset(key):
    return cam.offset(coord_mapping[key])

  def _scale(key):
    return cam.scale(coord_mapping[key])

  # dictionary
  rpc_dict = {
    "satId": "???",
    "bandId": "???",
    "SpecId": output_rational_order.name,
    "IMAGE": {
      # "errBias": ,
      # "errRand": ,
      "lineOffset": _offset("line"),
      "sampOffset": _offset("samp"),
      "latOffset": _offset("lat"),
      "longOffset": _offset("long"),
      "heightOffset": _offset("height"),
      "lineScale": _scale("line"),
      "sampScale": _scale("samp"),
      "latScale": _scale("lat"),
      "longScale": _scale("long"),
      "heightScale": _scale("height"),
      "lineNumCoef": _coefficients("lineNumCoef"),
      "lineDenCoef": _coefficients("lineDenCoef"),
      "sampNumCoef": _coefficients("sampNumCoef"),
      "sampDenCoef": _coefficients("sampDenCoef"),
    }
  }

  return rpc_dict

# add `as_dict` method to rational_camera
setattr(_vpgl.rational_camera, 'as_dict', rational_camera_as_dict)


def correct_rational_camera(cam, gt_offset_u, gt_offset_v, verbose=False):
  """
  Adds (delta-u, delta-v) offset to a rational camera.

  Parameters
  ----------
  cam : vpgl_rational_camera
    The camera to correct
  gt_offset_u : float
    The u offset
  gt_offset_v : float
    The v offset
  verbose : bool (optional)
    Tells VXL to print informative messages

  Returns
  -------
  New corrected vpgl_rational_camera

  Raises
  ------
  ValueError
    If cam is not a valid rational camera.
  """

  from ._vpgl import _correct_local_rational_camera, _correct_rational_camera

  if cam.type_name() == 'vpgl_local_rational_camera':
    return _correct_local_rational_camera(cam, gt_offset_u, gt_offset_v, verbose)
  elif cam.type_name() == 'vpgl_rational_camera':
    return _correct_rational_camera(cam, gt_offset_u, gt_offset_v, verbose)
  else:
    raise ValueError("Camera of type <{}> not a rational camera".format(cam.type_name()))


def _rsm_metadata_validate_dict(dct):
  '''
  This function searches through an rsm_metadata dictionary containing values
  protected by ``*_valid`` properties, sets any invalid properties to ``None``,
  and removed the ``*_valid`` properties.

  For example, the following rsm_metadata dictionary:

  .. code-block:: python

      dct = {
        "platform_name_valid": True,
        "platform_name": "A",
        "image_type_valid": False,
        "image_type": "garbage",
        ...
      }

  Would be transformed to:

  .. code-block:: python

      dct = {
        "platform_name": "A",
        "image_type": None,
        ...
      }

  '''

  # identify keys ending in `*_valid` (except for "any_valid")
  keys_valid = [key for key in dct.keys()
                if key.endswith('_valid') and key != "any_valid"]

  # keys requiring special treatment
  keys_invalidate = {
    "corners_valid": ["upper_left", "upper_right",
                      "lower_left", "lower_right",
                      "bounding_box", "footprint"],
    "image_corners_valid": ["min_image_corner", "max_image_corner"],
    "sun_azimuth_valid": ["sun_azimuth", "sun_azimuth_radians"],
    "sun_elevation_valid": ["sun_elevation", "sun_elevation_radians"],
  }

  # process keys
  for key_valid in keys_valid:
    is_valid = dct.pop(key_valid)
    if is_valid:
      continue

    keys = (keys_invalidate.get(key_valid) or
            [key_valid.replace('_valid', '')])
    for key in keys:
      dct[key] = None

  return dct

# add ``as_dict_validated`` method to rsm_metadata object
rsm_metadata.as_dict_validated = lambda x: _rsm_metadata_validate_dict(x.as_dict())
