import copy
import json
import numpy as np
import pickle
import unittest
import utils

from vxl import vgl
from vxl import vpgl
from vxl.contrib import acal


def json_serializer(obj):
  try:
    return str(obj)
  except err:
    raise TypeError("Type {} not serializable".format(type(obj))) from err


class AcalBase(utils.VxlBase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.default_data = dict()
    self.init_data = dict()

  def _data_after_init(self):
    data = copy.deepcopy(self.default_data)
    data = utils.update_nested_dict(data, getattr(self, 'init_data', {}))
    return data

  def _cls_instance(self, *args, **kwargs):
    # override to use different class creation method
    return self.cls(*args, **kwargs)

  def test_create(self):
    instance = self._cls_instance()
    self.assertIsInstance(instance, self.cls)
    self.assertAttributes(instance, self.default_data)

  @utils.skipUnlessAttr('init_data')
  def test_init(self):
    instance = self._cls_instance(**self.init_data)
    new_data = self._data_after_init()
    self.assertIsInstance(instance, self.cls)
    self.assertAttributes(instance, new_data)

  @utils.skipUnlessClassAttr('set')
  def test_set(self):
    instance = self._cls_instance()
    instance.set(**self.init_data)
    new_data = self._data_after_init()
    self.assertAttributes(instance, new_data)

  def test_equal(self):
    instance_A = self._cls_instance(**self.init_data)
    instance_B = self._cls_instance(**self.init_data)
    self.assertEqual(instance_A, instance_B)

  @utils.skipUnlessClassAttr('__getstate__')
  def test_pickle(self):
    instance_A = self._cls_instance(**self.init_data)
    instance_B = pickle.loads(pickle.dumps(instance_A))
    self.assertEqual(instance_A, instance_B)


class acal_f_params(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.f_params
    self.default_data = {
        'epi_dist_mul': 2.5,
        'max_epi_dist': 5.0,
        'F_similar_abcd_tol': 0.01,
        'F_similar_e_tol': 1.0,
        'ray_uncertainty_tol': 50.0,
        'min_num_matches': 5,
    }
    self.init_data = {
        'epi_dist_mul': 3.0,
        'max_epi_dist': 5.5,
        'F_similar_abcd_tol': 0.51,
        'F_similar_e_tol': 1.5,
        'ray_uncertainty_tol': 50.5,
        'min_num_matches': 6,
    }


class acal_corr(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.corr
    self.default_data = {
        'id': np.uint(-1),
        'pt': vgl.point_2d(-1, -1)
    }
    self.init_data = {
        'id': 10,
        'pt': vgl.point_2d(10.0, 20.0)
    }


class acal_match_pair(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_pair
    self.default_data = {
        'corr1': {'id': np.uint(-1), 'pt': vgl.point_2d(-1, -1)},
        'corr2': {'id': np.uint(-1), 'pt': vgl.point_2d(-1, -1)},
    }
    self.init_data = {
        'corr1': {'id': 10, 'pt': vgl.point_2d(10.0, 20.0)},
        'corr2': {'id': 15, 'pt': vgl.point_2d(15.0, 25.0)},
    }


class acal_match_params(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_params
    self.default_data = {
        'min_n_tracks': 3,
        'min_n_cams': 3,
        'max_proj_error': 1.0,
        'max_uncal_proj_error': 20.0,
    }
    self.init_data = {
        'min_n_tracks': 5,
        'min_n_cams': 10,
        'max_proj_error': 10.0,
        'max_uncal_proj_error': 30.0,
    }


class acal_match_node(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_node
    self.default_data = {
      'cam_id': 0,
      'children': [],
      'self_to_child_matches': [],
    }
    self.init_data = {
      'cam_id': 100,
    }


class acal_match_tree(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_tree
    self.default_data = {
      'root_id': 0,
      'min_n_tracks': 1,
    }
    self.init_data = {
      'root_id': 100,
    }


  def _construct_example(self):
    '''
    Example tree
    '''

    def _make_match_pair(id0, x0, y0, id1, x1, y1):
      return acal.match_pair(acal.corr(id0, vgl.point_2d(x0, y0)),
                             acal.corr(id1, vgl.point_2d(x1, y1)))

    # 12 -> 21 match pairs
    pairs_b = [
        ( 818, 188.987, 227.430,  1617, 278.163, 315.765),
        ( 983, 109.553, 284.380,  1834, 150.748, 361.344),
        (1019, 334.196, 286.747,  2022, 418.509, 393.688),
        (1073, 327.695, 300.218,  2163, 420.439, 419.478),
        (1075, 323.612, 301.511,  2159, 415.352, 420.685),
        (2395, 224.962, 100.704,   541, 282.362, 108.280),
    ]
    pairs_b = [_make_match_pair(*item) for item in pairs_b]

    # 21 -> 22 match pairs
    pairs_a = [
        (1988, 38.454, 388.516,  1378, 45.002, 400.503),
        (2047, 24.258, 406.217,  1468, 30.762, 419.905),
        (2159, 15.352, 420.685,  1427, 60.936, 405.096),
        (2163, 20.439, 419.478,  1429, 66.727, 403.056),
        (2207, 4.9441, 445.987,  1629, 3.7094, 476.475),
        (2308, 3.7507, 450.587,  1696, 2.5715, 481.109),
    ]
    pairs_a = [_make_match_pair(*item) for item in pairs_a]

    # tree structure (parent, child, pairs)
    root_id = 12
    tree_structure = [
        ( 12,  21, pairs_b),
        ( 21,  22, pairs_a),
        ( 12, 210, pairs_b),
        ( 21, 220, pairs_a),
        (210, 221, pairs_a),
        ( 21, 222, pairs_a),
    ]

    # create tree
    tree = acal.match_tree(root_id)

    # test root object
    root = tree.root
    self.assertEqual(root.cam_id, root_id,
                     "Empty tree, root node, incorrect id")
    self.assertEqual(root.is_root(), True,
                     "Empty tree, root node, incorrect is_root")
    self.assertEqual(root.is_leaf(), True,
                     "Empty tree, root node, incorrect is_leaf")

    # add single child
    child_data = tree_structure[0]
    tree.add_child_node(*child_data)

    self.assertEqual(len(root), 1,
                     "Single child, root node, incorrect size")
    self.assertEqual(root.is_root(), True,
                     "Single child, root node, incorrect is_root")
    self.assertEqual(root.is_leaf(), False,
                     "Single child, root node, incorrect is_leaf")

    child = root.children[0]
    self.assertEqual(child.parent_id(), child_data[0],
                     "Single child, child node, incorrect parent id")
    self.assertEqual(child.cam_id, child_data[1],
                     "Single child, child node, incorrect id")
    self.assertEqual(child.is_root(), False,
                     "Single child, child node, incorrect is_root")
    self.assertEqual(child.is_leaf(), True,
                     "Single child, child node, incorrect is_leaf")

    # add more children
    for item in tree_structure[1:]:
      tree.add_child_node(*item)

    root_children_ids = [item[1] for item in tree_structure
                         if item[0] == root_id]
    self.assertEqual(root.children_ids(), root_children_ids,
                     "Full tree, root node, incorrect children_ids")

    cam_ids = [i for item in tree_structure for i in item[0:2]]
    cam_ids = sorted(list(set(cam_ids)))
    self.assertEqual(tree.cam_ids(), cam_ids,
                     "Full tree, incorrect cam_ids")

    return tree


  def test_run(self):
    self._construct_example()


  def test_equal_after_run(self):
    tree_A = self._construct_example()
    tree_B = self._construct_example()
    self.assertEqual(tree_A, tree_B)


class acal_match_vertex(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_vertex
    self.default_data = {
      'cam_id': 0,
      'mark': False,
    }
    self.init_data = {
      'cam_id': 100,
    }


class acal_match_edge(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_edge
    self.default_data = {
      'id': 0,
      'v0_id': 0,
      'v1_id': 1,
      'matches': [],
    }
    self.init_data = {
      'id': 100,
      'v0_id': 100,
      'v1_id': 200,
    }

  def _cls_instance(self, **kwargs):

    if 'v0' not in kwargs:
      v0_id = kwargs.pop('v0_id', self.default_data['v0_id'])
      kwargs['v0'] = acal.match_vertex(cam_id = v0_id)

    if 'v1' not in kwargs:
      v1_id = kwargs.pop('v1_id', self.default_data['v1_id'])
      kwargs['v1'] = acal.match_vertex(cam_id = v1_id)

    if 'matches' not in kwargs:
      kwargs['matches'] = []

    return self.cls(**kwargs)


class acal_match_graph(AcalBase, unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.cls = acal.match_graph


  def _construct_example(self):
    '''
    Example graph inputs for two connected components
       component "A" = 4 images/cameras (index 0-3) with 4 correspondences
       component "B" = 3 images/cameras (index 4-6) with 3 correspondences
    For each image/camera, we project 3d points into the 2d image space
    to serve as feature correspondences.
    '''

    # helper: affine camera with more defaults
    def make_affine_camera(
        rayx, rayy, rayz,   # 3D ray direction
        upx = 0.0, upy = 0.0, upz = 1.0,  # 3D up direction
        ptx = 0.0, pty = 0.0, ptz = 0.0,  # 3D stare point
        u0 = 0, v0 = 0,  # stare point image projection
        su = 1.0, sv = 1.0,  # scaling
      ):
      return vpgl.affine_camera(vgl.vector_3d(rayx, rayy, rayz),
                                vgl.vector_3d(upx, upy, upz),
                                vgl.point_3d(ptx, pty, ptz),
                                u0, v0, su, sv)

    # helper: incidence matrix with points & pairs of cams
    def make_incidence(pts, pairs, cams, cam_offset = 0):

      incidence_list = []
      for i, j in pairs:
        vec = [acal.match_pair(acal.corr(k, cams[i].project(pt)),
                               acal.corr(k, cams[j].project(pt)))
               for k, pt in enumerate(pts)]
        incidence_list.append((i + cam_offset, j + cam_offset, vec))

      return incidence_list

    # component A
    rays = [(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)]
    camsA = [make_affine_camera(*r) for r in rays]

    pts = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1)]
    ptsA = [vgl.point_3d(*p) for p in pts]

    pairsA = [(0,1),(1,2),(2,3),(3,0)]

    cam_offset = 0
    incidenceA = make_incidence(ptsA, pairsA, camsA)

    # component B
    rays = [(1,1,0), (-1,1,0), (-1,-1,0)]
    camsB = [make_affine_camera(*r) for r in rays]

    pts = [(2,0,0), (0,2,0), (0,0,2)]
    ptsB = [vgl.point_3d(*p) for p in pts]

    pairsB = [(0,1),(1,2),(2,0)]

    cam_offset += len(camsA)
    incidenceB = make_incidence(ptsB, pairsB, camsB, cam_offset)

    # total system
    pts = ptsA + ptsB
    cams = {i: c for i, c in enumerate(camsA + camsB)}
    image_paths = {i: 'image{}.tif' for i in range(len(cams))}

    incidence_matrix = incidenceA + incidenceB
    incidence_matrix = {item[0]: {item[1]: item[2]} for item in incidence_matrix}
    # print(json.dumps(incidence_matrix, indent = 2, default = json_serializer))

    # create match graph
    match_graph = self._cls_instance()
    match_graph.acams = cams
    match_graph.image_paths = image_paths
    success = match_graph.load_incidence_matrix(incidence_matrix)
    self.assertTrue(success, "incidence matrix failed to load")

    # connected components
    match_graph.find_connected_components()
    components = match_graph.connected_components
    self.assertEqual(len(components), 2,
                     "incorrect number of connected components")
    self.assertEqual(len(components[0]), 4,
                     "incorrect size of connected component[0]")
    self.assertEqual(len(components[1]), 3,
                     "incorrect size of connected component[1]")

    # focus tracks
    match_graph.compute_focus_tracks()
    tracks = match_graph.focus_tracks
    self.assertEqual(len(tracks[0][0]), 4,
                     "incorrect size of focus_track[0][0]")
    self.assertEqual(len(tracks[1][4]), 3,
                     "incorrect size of focus_track[1][4]")

    # additional processing
    match_graph.compute_match_trees()
    match_graph.validate_match_trees_and_set_metric()

    return match_graph


  def test_run(self):
    self._construct_example()


  def test_equal_after_run(self):
    match_graph_A = self._construct_example()
    match_graph_B = self._construct_example()
    self.assertEqual(match_graph_A, match_graph_B)


if __name__ == '__main__':
  unittest.main()
