from ._bvxm_algo import *


def create_scene(scene_xml):
  """
  Create a bvxm_voxel_world from an xml file with the parameters.
  """

  # Local imports
  from ._bvxm_algo import _create_voxel_world
  from xml.etree import ElementTree

  # Read XML file
  scene_options = ElementTree.parse(scene_xml).getroot()

  # Parse parameters
  # TODO: Check that find() doesn't return None
  input_dir = scene_options.find('input_directory').get("value")
  lvcs_path = scene_options.find('lvcs').get("value")
  corner_x = float(scene_options.find('corner_x').get("value", 0.0))
  corner_y = float(scene_options.find('corner_y').get("value", 0.0))
  corner_z = float(scene_options.find('corner_z').get("value", 0.0))
  dimx = int(scene_options.find('voxel_dim_x').get("value", 10))
  dimy = int(scene_options.find('voxel_dim_y').get("value", 10))
  dimz = int(scene_options.find('voxel_dim_z').get("value", 10))
  vox_len = float(scene_options.find('voxel_length').get("value", 1))
  min_ocp_prob = float(scene_options.find('min_ocp_prob').get("value", 1e-5))
  max_ocp_prob = float(scene_options.find('max_ocp_prob').get("value", 1e-5))
  max_scale = int(scene_options.find('max_scale').get("value", 1))

  # Create voxel world, return it
  return _create_voxel_world(input_dir, lvcs_path,
                             corner_x, corner_y, corner_z,
                             dimx, dimy, dimz, vox_len,
                             min_ocp_prob, max_ocp_prob,
                             max_scale)
