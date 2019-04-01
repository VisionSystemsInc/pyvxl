#include "pyvpgl_file_formats.h"

#include <vpgl/file_formats/vpgl_geo_camera.h>
#include <vil/vil_load.h>
#include <vil/vil_image_resource.h>
#include <vil/vil_image_resource_sptr.h>

#include "../../pyvxl_util.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>
#include <vector>
#include <array>

namespace py = pybind11;

namespace pyvxl { namespace vpgl { namespace file_formats {

void wrap_vpgl_file_formats(py::module &m)
{
  // Geo- Camera definitions
  py::class_<vpgl_geo_camera>(m, "geo_camera")
    // Default methods
    .def(py::init<>())
    .def("__str__", streamToString<vpgl_geo_camera >)
    // Convert pixel coords (u,v) to a lon/lat pair
    .def("img_to_global",
      [](vpgl_geo_camera &G, double const u, double const v)
      {
        double lon, lat;
        G.img_to_global(u, v, lon, lat);
        return std::make_tuple(lon, lat);
      }
    );

  // Init from a Geotiff filename
  m.def("read_geo_camera",
    [](std::string filename)
    {
      vpgl_geo_camera* cam = new vpgl_geo_camera;
      vil_image_resource_sptr img = vil_load_image_resource(filename.c_str());
      vpgl_geo_camera::init_geo_camera(img, cam);
      return cam;
    },
    "A function to read a geo camera from a geotiff header."
  );
}
}}}


PYBIND11_MODULE(_vpgl_file_formats, m)
{
  m.doc() =  "Python bindings for the VIL file_formats computer vision libraries";

  pyvxl::vpgl::file_formats::wrap_vpgl_file_formats(m);
}
