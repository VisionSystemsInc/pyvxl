#include "pyvil.h"
#include <tuple>
#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace pyvxl { namespace vil {

template <class T>
vil_image_view<T>* image_from_buffer(py::array_t<T> b)
{
  py::buffer_info info = b.request();
  if (info.format != py::format_descriptor<T>::format()) {
    throw std::runtime_error("Incompatible scalar type");
  }
  if (info.ndim < 2) {
    throw std::runtime_error("Expecting at least 2 dimension");
  }
  if (info.ndim > 3) {
    throw std::runtime_error("Expecting at most 3 dimensions");
  }
  const size_t num_rows = info.shape[0];
  const size_t num_cols = info.shape[1];
  size_t num_planes = 1;
  if (info.ndim == 3) {
    num_planes = info.shape[2];
  }

  // in-place constructor
  vil_image_view<T> *img = new vil_image_view<T>(num_cols, num_rows, num_planes);
  size_t row_stride = info.strides[0]/sizeof(T);
  size_t col_stride = info.strides[1]/sizeof(T);
  size_t plane_stride = 0;
  if (info.ndim == 3) {
    plane_stride = info.strides[2]/sizeof(T);
  }
  for (size_t r=0; r<num_rows; ++r) {
    for (size_t c=0; c<num_cols; ++c) {
      for (size_t p=0; p<num_planes; ++p) {
        (*img)(c,r,p) = *(static_cast<T*>(info.ptr) + r*row_stride + c*col_stride + p*plane_stride);
      }
    }
  }

  return img;
}


template<class T>
T image_getitem(vil_image_view<T> const& img, std::tuple<size_t, size_t, size_t> pos)
{
  // TODO: wrap around
  size_t img_y = std::get<0>(pos);
  size_t img_x = std::get<1>(pos);
  size_t img_p = std::get<2>(pos);
  return img(img_x, img_y, img_p);
}

template<class T>
long image_len(vil_image_view<T> const& img)
{
  return img.ni()*img.nj()*img.nplanes();
}


template<class T>
std::tuple<size_t, size_t, size_t> image_view_shape(vil_image_view<T> const& img)
{
  return std::make_tuple(img.nj(), img.ni(), img.nplanes());
}

template<class T>
py::buffer_info get_image_buffer(vil_image_view<T> &img)
{
  int ndim = 2;
  std::vector<size_t> img_shape {img.nj(), img.ni()};
  std::vector<size_t> img_stride {static_cast<size_t>(img.jstep()),
                                  static_cast<size_t>(img.istep())};
  if (img.nplanes() > 1) {
    ndim = 3;
    img_shape.push_back(img.nplanes());
    img_stride.push_back(img.planestep());
  }
  return py::buffer_info(img.top_left_ptr(), sizeof(T),
                         py::format_descriptor<T>::format(),
                         ndim, img_shape,img_stride);
}

template<class T>
void wrap_vil_image_view(py::module &m, std::string const& class_name)
{
  py::class_<vil_image_view<T>, vil_image_view_base > (m, class_name.c_str(), py::buffer_protocol())
    .def(py::init<>())
    .def(py::init<unsigned int,unsigned int>())
    .def(py::init<unsigned int,unsigned int,unsigned int>())
    .def(py::init(&image_from_buffer<T>))
    .def("__len__", image_len<T>)
    .def("__getitem__", image_getitem<T>)
    .def_property_readonly("shape", &image_view_shape<T>)
    .def_buffer(get_image_buffer<T>);
}

template<class T>
vil_image_view<T> vil_load_wrapper(std::string const& filename)
{
  return vil_load(filename.c_str());
}

template<class T>
void vil_save_wrapper(vil_image_view<T> const& img, std::string const& filename)
{
  vil_save(img, filename.c_str());
}


void wrap_vil(py::module &m)
{
  // Need to provide a python wrapping for the base class so that 
  // the derived classes work correctly.  
  py::class_<vil_image_view_base>(m, "image_view_base");

  wrap_vil_image_view<unsigned char>(m, "image_view_byte");
  wrap_vil_image_view<float>(m, "image_view_float");
  // TODO: overload these so that they work on any pixel type
  m.def("load", &vil_load_wrapper<unsigned char>);
  m.def("save", &vil_save_wrapper<unsigned char>);
  // Lambda version of the above, in case that helps with the todo
  // m.def("load", [](std::string const& filename)
  // {
  //   return static_cast<vil_image_view<unsigned char>>(vil_load(filename.c_str()));
  // });
  // m.def("save", [](vil_image_view<unsigned char> img, std::string const& filename)
  // {
  //   return vil_save(img, filename.c_str());
  // });
}
}}
