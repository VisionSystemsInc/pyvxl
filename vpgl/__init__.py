from ._vpgl import *
from . import algo


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
