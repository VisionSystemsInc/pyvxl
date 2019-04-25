#include "pybrad.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <string>
#include <memory>
#include <sstream>
#include <vector>
#include <array>


#include <vil/vil_convert.h>
#include <vil/vil_image_view.h>
#include <vnl/vnl_math.h>
#include <brad/brad_image_atmospherics_est.h>
#include <brad/brad_image_metadata.h>
#include <brad/brad_calibration.h>

namespace py = pybind11;

namespace pyvxl { namespace brad {

vil_image_view<float> estimate_reflectance(vil_image_view<float> const& radiance_img,
                                           brad_image_metadata const& mdata,
                                           float mean_reflectance,
                                           bool average_airlight,
                                           bool is_normalize)
{

  if (radiance_img.pixel_format() != VIL_PIXEL_FORMAT_FLOAT) {
    throw std::invalid_argument("ERROR: vxl.contrib.brad.estimate_reflectance: expecting floating point radiance image\n");
  }

  unsigned int ni = radiance_img.ni();
  unsigned int nj = radiance_img.nj();
  unsigned int np = radiance_img.nplanes();
  if (mean_reflectance <= 0.0)
    is_normalize = false;

  auto* reflectance_img = new vil_image_view<float>(ni, nj, np);
  bool success = brad_estimate_reflectance_image(radiance_img, mdata, mean_reflectance, *reflectance_img, average_airlight, is_normalize);

  if (!success)
    throw std::runtime_error(std::string("ERROR: vxl.contrib.brad.estimate_reflectance: brad_estimate_reflectance_image failed.\n"));

  return *reflectance_img;
}


void wrap_brad(py::module &m)
{
  m.def("estimate_reflectance", &estimate_reflectance,
        py::arg("radiance"), py::arg("mdata"), py::arg("mean_reflectance"),
        py::arg("average_airlight"), py::arg("is_normalize"));

  m.def("radiometrically_calibrate", &brad_nitf_abs_radiometric_calibrate,
        py::arg("image"), py::arg("meta"));

  py::class_<image_time>(m, "image_time")
    .def_readwrite("year", &image_time::year)
    .def_readwrite("month", &image_time::month)
    .def_readwrite("day", &image_time::day)
    .def_readwrite("hour", &image_time::hour)
    .def_readwrite("min", &image_time::min)
    .def_readwrite("sec", &image_time::sec);

  py::class_<brad_image_metadata> (m, "brad_image_metadata")
    .def(py::init<>())
    .def(py::init<const std::string,std::string>())
    .def("parse", &brad_image_metadata::parse)
    .def("parse_from_meta_file", &brad_image_metadata::parse_from_meta_file)
    .def("same_time", &brad_image_metadata::same_time)
    .def("same_day", &brad_image_metadata::same_day)
    .def("time_minute_diff", &brad_image_metadata::time_minute_dif)
    .def("same_extent", &brad_image_metadata::same_extent)
    .def("read_band_dependent_gain_offset", &brad_image_metadata::read_band_dependent_gain_offset)
    .def("read_band_dependent_solar_irradiance", &brad_image_metadata::read_band_dependent_solar_irradiance)
    .def_readwrite("view_elevation", &brad_image_metadata::view_elevation_)
    .def_readwrite("view_azimuth", &brad_image_metadata::view_azimuth_)
    .def_readwrite("sun_elevation", &brad_image_metadata::sun_elevation_)
    .def_readwrite("sun_azimuth", &brad_image_metadata::sun_azimuth_)
    .def_readwrite("absolute_calibration", &brad_image_metadata::abscal_)
    .def_readwrite("effective_bandwidth", &brad_image_metadata::effect_band_width_)
    .def_readwrite("gains", &brad_image_metadata::gains_)
    .def_readwrite("offsets", &brad_image_metadata::offsets_)
    .def_readwrite("normal_sun_irradiance_values", &brad_image_metadata::normal_sun_irradiance_values_)
    .def_readwrite("sun_irradiance", &brad_image_metadata::sun_irradiance_)
    .def_readwrite("time", &brad_image_metadata::t_)
    .def_readwrite("number_of_bits", &brad_image_metadata::number_of_bits_)
    .def_readwrite("satellite_name", &brad_image_metadata::satellite_name_)
    .def_readwrite("cloud_coverage_percentage", &brad_image_metadata::cloud_coverage_percentage_)
    .def_readwrite("footprint", &brad_image_metadata::footprint_)
    .def_readwrite("lower_left", &brad_image_metadata::lower_left_)
    .def_readwrite("upper_right", &brad_image_metadata::upper_right_)
    .def_readwrite("camera_offset", &brad_image_metadata::cam_offset_)
    .def_readwrite("band", &brad_image_metadata::band_)
    .def_readwrite("number_bands", &brad_image_metadata::n_bands_)
    .def_readwrite("gsd", &brad_image_metadata::gsd_)
    ;
}

}}

PYBIND11_MODULE(_brad, m)
{
  m.doc() =  "Python bindings for the VXL BRAD computer vision library";

  pyvxl::brad::wrap_brad(m);

}

