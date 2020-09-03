#include "pybvxm.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <bvxm/bvxm_edge_util.h>
#include <bvxm/bvxm_voxel_world.h>
#include <bvxm/bvxm_world_params.h>
#include <vpgl/vpgl_utm.h>

namespace py = pybind11;

namespace pyvxl { namespace bvxm {


// generates an ortho camera from the scene bounding box, GSD of the image is 1 meter
vpgl_geo_camera create_ortho_camera(bvxm_voxel_world const& world, bool is_utm = false)
{
  // generate vpgl_geo_camera for the scene
  bvxm_world_params_sptr params = world.get_params();
  vgl_box_3d<double> box = params->world_box_local();
  vgl_point_3d<float> corner = params->corner();
  vgl_point_3d<float> upper_left(corner.x(), (float)(corner.y() + box.height()), corner.z());
  vgl_point_3d<float> lower_right((float)(corner.x()+box.width()), corner.y(), corner.z());
  float voxel_length = params->voxel_length();

  vpgl_lvcs_sptr lvcs = params->lvcs();
  double lat, lon, elev;
  lvcs->get_origin(lat, lon, elev);
  std::cout << " lvcs origin: " << lat << " " << lon << " " << elev << std::endl;

  // determine the upper left corner to use a vpgl_geo_cam, WARNING: assumes that the world is compass-alinged
  double upper_left_lon, upper_left_lat, upper_left_elev;
  lvcs->local_to_global(upper_left.x(), upper_left.y(), upper_left.z(), vpgl_lvcs::wgs84, upper_left_lon, upper_left_lat, upper_left_elev);
  std::cout << "upper left corner in the image is: " << upper_left_lon << " lat: " << upper_left_lat << std::endl;

  double lower_right_lon, lower_right_lat, lower_right_elev;
  lvcs->local_to_global(lower_right.x(), lower_right.y(), lower_right.z(), vpgl_lvcs::wgs84, lower_right_lon, lower_right_lat, lower_right_elev);
  std::cout << "lower right corner in the image is: " << lower_right_lon << " lat: " << lower_right_lat << std::endl;

  vnl_matrix<double> trans_matrix(4,4,0.0);

  if (is_utm) {
    double scale_x = voxel_length;
    double scale_y = -1*voxel_length;
    // transfer upper left corner to utm system
    vpgl_utm utm;
    double upper_left_x, upper_left_y;
    int utm_zone;
    utm.transform(upper_left_lat, upper_left_lon, upper_left_x, upper_left_y, utm_zone);
    std::cout << "upper left in utm = " << upper_left_x << " x " << upper_left_y << std::endl;
    std::cout << "scale_x = " << scale_x << " scale_y = " << scale_y << std::endl;
    trans_matrix[0][0] = scale_x;
    trans_matrix[1][1] = scale_y;
    trans_matrix[0][3] = upper_left_x;
    trans_matrix[1][3] = upper_left_y;
    vpgl_geo_camera cam(trans_matrix, lvcs);
    unsigned northing = 1;
    if (upper_left_lat < 0 && lower_right_lat < 0)
      northing = 0;
    if (upper_left_lat*lower_right_lat < 0)
      std::cout << "warning: scene world crosses the Equator" << std::endl;
    cam.set_utm(utm_zone,northing);
    cam.set_scale_format(true);
    return cam;
  }
  else {
    auto ni = (unsigned)std::ceil(box.width() / voxel_length);
    auto nj = (unsigned)std::ceil(box.height()/ voxel_length);
    //trans_matrix[0][0] = (lower_right_lon-lon)/ni; trans_matrix[1][1] = -(upper_left_lat-lat)/nj;
    // lvcs origin is not necessarily one of the corners of the scene
    trans_matrix[0][0] = (lower_right_lon-upper_left_lon)/ni; trans_matrix[1][1] = -(upper_left_lat-lower_right_lat)/nj;
    trans_matrix[0][3] = upper_left_lon; trans_matrix[1][3] = upper_left_lat;
    vpgl_geo_camera cam(trans_matrix, lvcs);
    cam.set_scale_format(true);
    return cam;
  }


}


void wrap_bvxm(py::module &m)
{
  py::class_<bvxm_voxel_world> (m, "voxel_world")
    .def(py::init<>())
    .def("get_box",
      [](bvxm_voxel_world const& world)
      {
        bvxm_world_params_sptr params = world.get_params();

        vpgl_lvcs_sptr lvcs = params->lvcs();

        double lower_left_lon, lower_left_lat, lower_left_elev, upper_right_lon, upper_right_lat, upper_right_elev;
        lvcs->local_to_global(params->corner().x(), params->corner().y(), params->corner().z(), vpgl_lvcs::wgs84,
                              lower_left_lon, lower_left_lat, lower_left_elev);
        double dimx = params->num_voxels().x() * params->voxel_length();
        double dimy = params->num_voxels().y() * params->voxel_length();
        double dimz = params->num_voxels().z() * params->voxel_length();
        lvcs->local_to_global(params->corner().x() + dimx, params->corner().y() + dimy, params->corner().z() + dimz, vpgl_lvcs::wgs84,
                              upper_right_lon, upper_right_lat, upper_right_elev);

        return std::make_tuple(lower_left_lon, lower_left_lat, lower_left_elev,
                              upper_right_lon, upper_right_lat, upper_right_elev);
      })

    .def("get_local_box",
      [](bvxm_voxel_world const& world)
      {
        bvxm_world_params_sptr params = world.get_params();

        double min_x, min_y, max_x, max_y, voxel_size, min_z, max_z;
        min_x = params->corner().x();
        min_y = params->corner().y();
        min_z = params->corner().z();
        max_x = params->num_voxels().x() * params->voxel_length() + params->corner().x();
        max_y = params->num_voxels().y() * params->voxel_length() + params->corner().y();
        max_z = params->num_voxels().z() * params->voxel_length() + params->corner().z();
        voxel_size = params->voxel_length();

        return std::make_tuple(min_x, min_y, max_x, max_y, voxel_size, min_z, max_z);
      })

    .def("get_origin",
      [](bvxm_voxel_world const& world)
      {
        bvxm_world_params_sptr params = world.get_params();

        vpgl_lvcs_sptr lvcs = params->lvcs();

        double lower_left_lon, lower_left_lat, gz;
        lvcs->local_to_global(params->corner().x(), params->corner().y(), params->corner().z(), vpgl_lvcs::wgs84, lower_left_lon, lower_left_lat, gz);

        return std::make_tuple(lower_left_lon, lower_left_lat, gz);
      })

    .def("get_model_dir",
      [](bvxm_voxel_world const& world)
      {
        bvxm_world_params_sptr params = world.get_params();

        std::string model_dir = params->model_dir();

        return model_dir;
      });


  m.def("create_ortho_camera", &create_ortho_camera,
        py::arg("scene"), py::arg("is_utm") = false);

  m.def("detect_edges", &bvxm_edge_util::detect_edges,
        py::arg("input_img"), py::arg("noise_multiplier"), py::arg("smooth"),
        py::arg("automatic_threshold"), py::arg("junctionp"), py::arg("aggressive_junction_closure"));
}

}}

PYBIND11_MODULE(_bvxm, m)
{
  m.doc() = "Python bindings for the VXL BVXM computer vision library";

  pyvxl::bvxm::wrap_bvxm(m);
}

