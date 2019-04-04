#include "pybrad.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>
#include <vector>
#include <array>

#include <vil/vil_image_view.h>
#include <brad/brad_image_atmospherics_est.h>
#include <brad/brad_image_metadata.h>

namespace py = pybind11;

namespace pyvxl { namespace brad {

void wrap_brad(py::module &m)
{
  m.def("estimate_reflectance", &brad_estimate_reflectance_image);

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

