#include "pybrip.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <sstream>
#include <stdexcept>

#include <brip/brip_vil_nitf_ops.h>
#include <vil/vil_image_view.h>

namespace py = pybind11;

namespace pyvxl { namespace brip {

vil_image_view<vxl_byte> _truncate_nitf_image_to_byte(vil_image_view<vxl_uint_16> input_img, bool is_scale)
{
  // truncate the input 16 bits image to a byte image by ignoring the most significant 5 bits and less significant 3 bits
  std::cout << "vxl.contrib.brip.truncate_nitf_image truncating to byte image" << std::endl;

  unsigned ni = input_img.ni();
  unsigned nj = input_img.nj();
  unsigned np = input_img.nplanes();
  vil_image_view<vxl_byte> out_img(ni, nj, np);

  if (is_scale) {
    if (!brip_vil_nitf_ops::scale_nitf_bits(input_img, out_img)) {
      std::ostringstream buffer;
      buffer << "vxl.contrib.brip.truncate_nitf_image_to_byte: scaling nitf image from " << input_img.pixel_format() << " to " << out_img.pixel_format() << " failed!" << std::endl;
      throw std::runtime_error(buffer.str());
    }
  }
  else {
    if (!brip_vil_nitf_ops::truncate_nitf_bits(input_img, out_img)) {
      std::ostringstream buffer;
      buffer << "vxl.contrib.brip.truncate_nitf_image_to_byte: truncating nitf image from " << input_img.pixel_format() << " to " << out_img.pixel_format() << " failed!" << std::endl;
      throw std::runtime_error(buffer.str());
    }
  }

  return out_img;

}

vil_image_view<vxl_uint_16> _truncate_nitf_image_to_short(vil_image_view<vxl_uint_16> input_img)
{
  std::cout << "vxl.contrib.brip.truncate_nitf_image truncating to short image" << std::endl;

  unsigned ni = input_img.ni();
  unsigned nj = input_img.nj();
  unsigned np = input_img.nplanes();
  vil_image_view<vxl_uint_16> out_img(ni, nj, np);

  // truncate the input image by ignoring the most significant 5 bits
  if (!brip_vil_nitf_ops::truncate_nitf_bits(input_img, out_img)) {
    std::ostringstream buffer;
    buffer << "vxl.contrib.brip.truncate_nitf_image_to_short: truncating nitf image from " << input_img.pixel_format() << " to " << out_img.pixel_format() << " failed!" << std::endl;
    throw std::runtime_error(buffer.str());
  }

  return out_img;
}


void wrap_brip(py::module &m)
{
  m.def("_truncate_nitf_image_to_byte", &_truncate_nitf_image_to_byte,
        py::call_guard<py::gil_scoped_release>(),
        py::arg("input_img"), py::arg("is_scale"));
  m.def("_truncate_nitf_image_to_short", &_truncate_nitf_image_to_short,
        py::call_guard<py::gil_scoped_release>(),
        py::arg("input_img"));
}

}}


PYBIND11_MODULE(_brip, m)
{
  m.doc() = "Python bindings for the VXL BRIP computer vision library";

  pyvxl::brip::wrap_brip(m);
}

