#include "pybvxm_algo.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <bvxm/bvxm_voxel_world.h>
#include <bvxm/bvxm_world_params.h>
#include <bvxm/algo/bvxm_create_scene_xml.h>
#include <vul/vul_file.h>

namespace py = pybind11;

namespace pyvxl { namespace bvxm { namespace algo {

bvxm_voxel_world create_voxel_world(std::string const& vox_dir, std::string const& lvcs_path,
                        float corner_x = 0.0f, float corner_y = 0.0f, float corner_z = 0.0f,
                        unsigned int dimx = 10, unsigned int dimy = 10, unsigned int dimz = 10,
                        float vox_len = 1.0f, float min_ocp_prob = 1e-5f, float max_ocp_prob = 1- 1e-5f,
                        unsigned int max_scale = 1)
{
  if (!vul_file::is_directory(vox_dir) || !vul_file::exists(vox_dir)) {
    std::ostringstream buffer;
    buffer << "In bvxm.create_voxel_world -- input directory "  << vox_dir << "is not valid!\n";
    throw std::invalid_argument(buffer.str());
  }

  vgl_point_3d<float> corner(corner_x, corner_y, corner_z);

  vgl_vector_3d<unsigned int> voxel_dims(dimx, dimy, dimz);

  vpgl_lvcs_sptr lvcs = new vpgl_lvcs();
  if (lvcs_path != "") {
    std::ifstream is(lvcs_path.c_str());
    if (!is)
    {
      std::ostringstream buffer;
      buffer << " Error opening file  " << lvcs_path << std::endl;
      throw std::invalid_argument(buffer.str());
    }
    lvcs->read(is);
  }

  bvxm_world_params_sptr params = new bvxm_world_params();
  params->set_params(vox_dir, corner, voxel_dims, vox_len, lvcs, min_ocp_prob, max_ocp_prob, max_scale);

  bvxm_voxel_world vox_world;
  vox_world.set_params(params);

  return vox_world;
}


void wrap_bvxm_algo(py::module &m)
{

  m.def("_create_voxel_world", &create_voxel_world,
        py::arg("vox_dir"), py::arg("lvcs_path"),
        py::arg("corner_x") = 0.0f, py::arg("corner_y") = 0.0f, py::arg("corner_z") = 0.0f,
        py::arg("dimx") = 10, py::arg("dimy") = 10, py::arg("dimz") = 10,
        py::arg("vox_len") = 1.0f, py::arg("min_ocp_prob") = 1e-5f, py::arg("max_ocp_prob") = 1- 1e-5f,
        py::arg("max_scale") = 1);

  /* m.def("create_scene", &bvxm_create_scene); */
  m.def("create_scene_large_scale", &bvxm_create_scene_xml_large_scale,
        py::arg("bbox"), py::arg("scene_root"),
        py::arg("world_dir"), py::arg("dem_folder"),
        py::arg("world_size_in")=500.0, py::arg("voxel_size")=1.0,
        py::arg("height_diff")=120.0, py::arg("height_sub")=25.0,
        py::arg("extension")=500.0, py::arg("land_folder") = "");
}

}}}

PYBIND11_MODULE(_bvxm_algo, m)
{
  m.doc() = "Python bindings for the VXL BVXM algo computer vision library";

  pyvxl::bvxm::algo::wrap_bvxm_algo(m);
}

