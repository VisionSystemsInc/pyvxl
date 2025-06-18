#include "pyvil.h"
#include <tuple>
#include <vil/vil_convert.h>
#include <vil/vil_crop.h>
#include <vil/vil_image_view.h>
#include <vil/vil_image_view_base.h>
#include <vil/vil_load.h>
#include <vil/vil_math.h>
#include <vil/vil_new.h>
#include <vil/vil_pixel_format.h>
#include <vil/vil_save.h>
#include <vil/vil_smart_ptr.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "pyvxl_holder_types.h"

namespace py = pybind11;

namespace pyvxl { namespace vil {

/* This is a "trampoline" helper class for the virtual vil_image_resource
 * class, which redirects virtual calls back to Python */
class PyImageResource : public vil_image_resource {
public:
  using vil_image_resource::vil_image_resource; /* Inherit constructors */

  unsigned nplanes() const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      unsigned,           /* Return type */
      vil_image_resource, /* Parent class */
      nplanes,            /* Name of function */
                          /* No arguments, the trailing comma
                            * is needed for some compilers */
    );
  }

  unsigned ni() const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      unsigned,           /* Return type */
      vil_image_resource, /* Parent class */
      ni,                 /* Name of function */
                          /* No arguments, the trailing comma
                           * is needed for some compilers */
    );
  }

  unsigned nj() const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      unsigned,           /* Return type */
      vil_image_resource, /* Parent class */
      nj,                 /* Name of function */
                          /* No arguments, the trailing comma
                            * is needed for some compilers */
    );
  }

  enum vil_pixel_format pixel_format() const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      enum vil_pixel_format,  /* Return type */
      vil_image_resource,     /* Parent class */
      pixel_format,           /* Name of function */
                              /* No arguments, the trailing comma
                              * is needed for some compilers */
    );
  }

  vil_image_view_base_sptr get_view(unsigned i0, unsigned n_i,
                                    unsigned j0, unsigned n_j) const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD(
      vil_image_view_base_sptr,  /* Return type */
      vil_image_resource,        /* Parent class */
      get_view,                  /* Name of function */
      i0, n_i, j0, n_j           /* Arguments */
    );
  }

  vil_image_view_base_sptr get_copy_view(unsigned i0, unsigned n_i,
                                          unsigned j0, unsigned n_j) const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      vil_image_view_base_sptr,  /* Return type */
      vil_image_resource,        /* Parent class */
      get_copy_view,                  /* Name of function */
      i0, n_i, j0, n_j           /* Arguments */
    );
  }

  bool put_view(const vil_image_view_base& im, unsigned i0, unsigned j0) override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      bool,                /* Return type */
      vil_image_resource,  /* Parent class */
      put_view,            /* Name of function */
      im, i0, j0           /* Arguments */
    );
  }

  bool put_view(const vil_image_view_base& im) override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD(
      bool,                /* Return type */
      vil_image_resource,  /* Parent class */
      put_view,            /* Name of function */
      im                   /* Arguments */
    );
  }

  bool view_fits(const vil_image_view_base& im, unsigned i0, unsigned j0) override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD(
      bool,                /* Return type */
      vil_image_resource,  /* Parent class */
      view_fits,           /* Name of function */
      im, i0, j0           /* Arguments */
    );
  }

  char const* file_format() const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD(
      char const*,         /* Return type */
      vil_image_resource,  /* Parent class */
      file_format,         /* Name of function */
                            /* No arguments, the trailing comma
                            * is needed for some compilers */
    );
  }

  bool get_property(char const* tag, void* property_value = nullptr) const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      bool,                /* Return type */
      vil_image_resource,  /* Parent class */
      get_property,        /* Name of function */
      tag, property_value  /* Arguments */
    );
  }

};


/* This is a "trampoline" helper class for the virtual vil_image_view_base
 * class, which redirects virtual calls back to Python */
class PyBaseView : public vil_image_view_base {
public:
  using vil_image_view_base::vil_image_view_base; /* Inherit constructors */

  void set_size(unsigned width, unsigned height) override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      void,                 /* Return type */
      vil_image_view_base,  /* Parent class */
      set_size,             /* Name of function */
      width, height         /* Arguments */
    );
  }

  void set_size(unsigned width, unsigned height, unsigned n_planes) override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      void,                    /* Return type */
      vil_image_view_base,     /* Parent class */
      set_size,                /* Name of function */
      width, height, n_planes  /* Arguments */
    );
  }

  void print(std::ostream& os) const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      void,                 /* Return type */
      vil_image_view_base,  /* Parent class */
      print,                /* Name of function */
      os                    /* Arguments */
    );
  }

  std::string is_a() const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      std::string,          /* Return type */
      vil_image_view_base,  /* Parent class */
      is_a,                 /* Name of function */
                            /* No arguments, the trailing comma
                             * is needed for some compilers */
    );
  }

  enum vil_pixel_format pixel_format() const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      enum vil_pixel_format,  /* Return type */
      vil_image_view_base,    /* Parent class */
      pixel_format,           /* Name of function */
                              /* No arguments, the trailing comma
                               * is needed for some compilers */
    );
  }

  bool is_class(std::string const& s) const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD(
      bool,                 /* Return type */
      vil_image_view_base,  /* Parent class */
      is_class,             /* Name of function */
      s                     /* Arguments */
    );
  }

};


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

std::pair<uint64_t, bool> __multiply_with_overflow(size_t a, size_t b) {
  uint64_t ret = a * b;
  bool no_overflow;
  if (a != 0 && (ret / a) != b) {
    no_overflow = false;
  } else {
    no_overflow = true;
  }
  return std::pair<uint64_t, bool>(ret, no_overflow);
}

template<class T>
uint64_t image_len(vil_image_view<T> const& img)
{
  std::pair<uint64_t, bool> single_channel_area = __multiply_with_overflow(img.ni(), img.nj());
  if (single_channel_area.second) {
    std::pair<uint64_t, bool> all_channels_area = __multiply_with_overflow(single_channel_area.first, img.nplanes());
    if (all_channels_area.second) {
      return all_channels_area.first;
    }
  }
  std::ostringstream buffer;
  buffer << "Image view length overflowed! Image view had dimensions:"
    << std::endl << "ni: " << img.ni() << std::endl << "nj: " << img.nj()
    << std::endl << "nplanes: " << img.nplanes() << std::endl;
  throw std::runtime_error(buffer.str());
}

uint64_t resource_len(vil_image_resource const& img)
{
  std::pair<uint64_t, bool> single_channel_area = __multiply_with_overflow(img.ni(), img.nj());
  if (single_channel_area.second) {
    std::pair<uint64_t, bool> all_channels_area = __multiply_with_overflow(single_channel_area.first, img.nplanes());
    if (all_channels_area.second) {
      return all_channels_area.first;
    }
  }
  std::ostringstream buffer;
  buffer << "Image resource length overflowed! Image resource had dimensions:"
    << std::endl << "ni: " << img.ni() << std::endl << "nj: " << img.nj()
    << std::endl << "nplanes: " << img.nplanes() << std::endl;
  throw std::runtime_error(buffer.str());
}


template<class T>
std::tuple<size_t, size_t, size_t> image_view_shape(vil_image_view<T> const& img)
{
  return std::make_tuple(img.nj(), img.ni(), img.nplanes());
}

std::tuple<size_t, size_t, size_t> image_resource_shape(vil_image_resource const& img)
{
  return std::make_tuple(img.nj(), img.ni(), img.nplanes());
}

template<class T>
py::buffer_info get_image_buffer(vil_image_view<T> &img)
{
  size_t nbytes = sizeof(T);
  int ndim = 2;
  std::vector<size_t> img_shape {img.nj(), img.ni()};
  std::vector<size_t> img_stride {static_cast<size_t>(img.jstep()) * nbytes,
                                  static_cast<size_t>(img.istep()) * nbytes};
  if (img.nplanes() > 1) {
    ndim = 3;
    img_shape.push_back(img.nplanes());
    img_stride.push_back(static_cast<size_t>(img.planestep()) * nbytes);
  }
  return py::buffer_info(img.top_left_ptr(), sizeof(T),
                         py::format_descriptor<T>::format(),
                         ndim, img_shape,img_stride);
}

template<class T>
void wrap_vil_image_view(py::module &m, std::string const& class_name)
{
  py::class_<vil_image_view<T>, vil_image_view_base /* <- Parent */ > (m, class_name.c_str(), py::buffer_protocol())
    .def(py::init<>())
    .def(py::init<unsigned int,unsigned int>())
    .def(py::init<unsigned int,unsigned int,unsigned int>())
    .def(py::init(&image_from_buffer<T>))
    .def("__len__", image_len<T>)
    .def("__getitem__", image_getitem<T>)
    .def_property_readonly("shape", &image_view_shape<T>)
    .def_buffer(get_image_buffer<T>)
    .def("deep_copy", &vil_image_view<T>::deep_copy);
}

vil_image_view<unsigned char> load_byte(std::string filename)
{
  auto view = vil_convert_cast(vxl_byte(), vil_load(filename.c_str()));
  if (!view)
    throw std::runtime_error("Failed to load byte image " + filename);
  return view;
}

vil_image_view<unsigned short int> load_short(std::string filename)
{
  auto view = vil_convert_cast(vxl_uint_16(), vil_load(filename.c_str()));
  if (!view)
    throw std::runtime_error("Failed to load short image " + filename);
  return view;
}

vil_image_view<float> load_float(std::string filename)
{
  auto view = vil_convert_cast(float(), vil_load(filename.c_str()));
  if (!view)
    throw std::runtime_error("Failed to load float image " + filename);
  return view;
}

vil_image_view<int> load_int(std::string filename)
{
  auto view = vil_convert_cast(int(), vil_load(filename.c_str()));
  if (!view)
    throw std::runtime_error("Failed to load int image " + filename);
  return view;
}

vil_image_resource_sptr vil_load_image_resource_wrapper(std::string const& filename)
{
  return vil_load_image_resource(filename.c_str());
}

template <class T>
void vil_save_wrapper(vil_image_view<T> const& img, std::string const& filename)
{
  bool result = vil_save(img, filename.c_str());

  if ( !result ) {
    std::ostringstream buffer;
    buffer << "Failed to save image to " << filename << std::endl;
    throw std::runtime_error(buffer.str());
  }
}

template <class T>
double vil_image_sum_wrapper(vil_image_view<T> const& image, unsigned p = 0)
{
  // compute the sum of elements in plane p of image

  // store result in reference
  double sum = 0.0;
  vil_math_sum(sum, image, p);

  return sum;
}

template <class T>
vil_image_view<unsigned char> vil_stretch_image_to_byte_wrapper(vil_image_view<T> const& image, float min_limit, float max_limit)
{
  // convert input image to float
  vil_image_view<float> fimage;

  // if the src imagery is already of type float then vil_convert_cast simply does a
  // shallow copy. so as not to modify the original imagery, make a deep copy
  if( vil_pixel_format_component_format(image.pixel_format()) ==
      VIL_PIXEL_FORMAT_FLOAT) {
    fimage.deep_copy(image);
  }
  else {
    vil_convert_cast(image, fimage);
  }

  unsigned ni = fimage.ni(), nj = fimage.nj(), np = fimage.nplanes();

  // stretch out float image to cover full [0,255] range
  float scale = 255.0f / (max_limit - min_limit);
  for(unsigned j = 0; j<nj; ++j)
    for(unsigned i = 0; i<ni; ++i)
      for(unsigned p = 0; p<np; ++p){
        float v = scale * (fimage(i,j,p) - min_limit);
        if(v > 255.0f) v = 255.0f; if(v < 0.0f) v = 0.0f;
        fimage(i,j,p) = v;
      }

  // convert to byte image
  vil_image_view<unsigned char> byte_image;
  vil_convert_cast(fimage, byte_image);

  return byte_image;
}

template <class T>
vil_image_view<unsigned short int> vil_stretch_image_to_short_wrapper(vil_image_view<T> const& image, float min_limit, float max_limit)
{
  // convert input image to float
  vil_image_view<float> fimage;

  // if the src imagery is already of type float then vil_convert_cast simply does a
  // shallow copy. so as not to modify the original imagery, make a deep copy
  if( vil_pixel_format_component_format(image.pixel_format()) ==
      VIL_PIXEL_FORMAT_FLOAT) {
    fimage.deep_copy(image);
  }
  else {
    vil_convert_cast(image, fimage);
  }

  unsigned ni = fimage.ni(), nj = fimage.nj(), np = fimage.nplanes();

  float scale = 65536.0f / (max_limit - min_limit);
  for(unsigned j = 0; j<nj; ++j)
    for(unsigned i = 0; i<ni; ++i)
      for(unsigned p = 0; p<np; ++p){
        float v = scale * (fimage(i,j,p) - min_limit);
        if(v > 65536.0f) v= 65536.0f; if(v < 0.0f) v = 0.0f;
        fimage(i,j,p) = v;
      }

  // convert to short image
  vil_image_view<unsigned short int> short_image;
  vil_convert_cast(fimage, short_image);

  return short_image;
}

template <class T>
vil_image_view<float> vil_stretch_image_to_float_wrapper(vil_image_view<T> const& image, float min_limit, float max_limit)
{

  // convert input image to float
  vil_image_view<float> fimage;

  // if the src imagery is already of type float then vil_convert_cast simply does a
  // shallow copy. so as not to modify the original imagery, make a deep copy
  if( vil_pixel_format_component_format(image.pixel_format()) ==
      VIL_PIXEL_FORMAT_FLOAT) {
    fimage.deep_copy(image);
  }
  else {
    vil_convert_cast(image, fimage);
  }

  unsigned ni = fimage.ni(), nj = fimage.nj(), np = fimage.nplanes();

  float scale = 1.0f / (max_limit - min_limit);
  for(unsigned j = 0; j<nj; ++j)
    for(unsigned i = 0; i<ni; ++i)
      for(unsigned p = 0; p<np; ++p){
        float v = scale * (fimage(i,j,p) - min_limit);
        if(v > 1.0f) v = 1.0f; if(v < 0.0f) v=0.0f;
        fimage(i,j,p) = v;
      }

  return fimage;
}


template <class T>
std::tuple<T, T> vil_image_range_wrapper(vil_image_view<T> const& img) {
  T min_val, max_val;
  vil_math_value_range(img, min_val, max_val);
  return std::make_tuple(min_val, max_val);
}


void wrap_vil(py::module &m)
{

  py::enum_<vil_pixel_format>(m, "pixel_format", py::arithmetic(), "The VXL pixel type")
    .value("VIL_PIXEL_FORMAT_UNKNOWN", vil_pixel_format::VIL_PIXEL_FORMAT_UNKNOWN)
#if VXL_HAS_INT_64
    .value("VIL_PIXEL_FORMAT_UINT_64", vil_pixel_format::VIL_PIXEL_FORMAT_UINT_64)
    .value("VIL_PIXEL_FORMAT_INT_64", vil_pixel_format::VIL_PIXEL_FORMAT_INT_64)
#endif
    .value("VIL_PIXEL_FORMAT_UINT_32", vil_pixel_format::VIL_PIXEL_FORMAT_UINT_32)
    .value("VIL_PIXEL_FORMAT_INT_32", vil_pixel_format::VIL_PIXEL_FORMAT_INT_32)
    .value("VIL_PIXEL_FORMAT_UINT_16", vil_pixel_format::VIL_PIXEL_FORMAT_UINT_16)
    .value("VIL_PIXEL_FORMAT_INT_16", vil_pixel_format::VIL_PIXEL_FORMAT_INT_16)
    .value("VIL_PIXEL_FORMAT_BYTE", vil_pixel_format::VIL_PIXEL_FORMAT_BYTE)
    .value("VIL_PIXEL_FORMAT_SBYTE", vil_pixel_format::VIL_PIXEL_FORMAT_SBYTE)
    .value("VIL_PIXEL_FORMAT_FLOAT", vil_pixel_format::VIL_PIXEL_FORMAT_FLOAT)
    .value("VIL_PIXEL_FORMAT_DOUBLE", vil_pixel_format::VIL_PIXEL_FORMAT_DOUBLE)
    .value("VIL_PIXEL_FORMAT_BOOL", vil_pixel_format::VIL_PIXEL_FORMAT_BOOL)
#if VXL_HAS_INT_64
    .value("VIL_PIXEL_FORMAT_RGB_UINT_64", vil_pixel_format::VIL_PIXEL_FORMAT_RGB_UINT_64)
    .value("VIL_PIXEL_FORMAT_RGB_INT_64", vil_pixel_format::VIL_PIXEL_FORMAT_RGB_INT_64)
#endif
    .value("VIL_PIXEL_FORMAT_RGB_UINT_32", vil_pixel_format::VIL_PIXEL_FORMAT_RGB_UINT_32)
    .value("VIL_PIXEL_FORMAT_RGB_INT_32", vil_pixel_format::VIL_PIXEL_FORMAT_RGB_INT_32)
    .value("VIL_PIXEL_FORMAT_RGB_UINT_16", vil_pixel_format::VIL_PIXEL_FORMAT_RGB_UINT_16)
    .value("VIL_PIXEL_FORMAT_RGB_INT_16", vil_pixel_format::VIL_PIXEL_FORMAT_RGB_INT_16)
    .value("VIL_PIXEL_FORMAT_RGB_BYTE", vil_pixel_format::VIL_PIXEL_FORMAT_RGB_BYTE)
    .value("VIL_PIXEL_FORMAT_RGB_SBYTE", vil_pixel_format::VIL_PIXEL_FORMAT_RGB_SBYTE)
    .value("VIL_PIXEL_FORMAT_RGB_FLOAT", vil_pixel_format::VIL_PIXEL_FORMAT_RGB_FLOAT)
    .value("VIL_PIXEL_FORMAT_RGB_DOUBLE", vil_pixel_format::VIL_PIXEL_FORMAT_RGB_DOUBLE)
#if VXL_HAS_INT_64
    .value("VIL_PIXEL_FORMAT_RGBA_UINT_64", vil_pixel_format::VIL_PIXEL_FORMAT_RGBA_UINT_64)
    .value("VIL_PIXEL_FORMAT_RGBA_INT_64", vil_pixel_format::VIL_PIXEL_FORMAT_RGBA_INT_64)
#endif
    .value("VIL_PIXEL_FORMAT_RGBA_UINT_32", vil_pixel_format::VIL_PIXEL_FORMAT_RGBA_UINT_32)
    .value("VIL_PIXEL_FORMAT_RGBA_INT_32", vil_pixel_format::VIL_PIXEL_FORMAT_RGBA_INT_32)
    .value("VIL_PIXEL_FORMAT_RGBA_UINT_16", vil_pixel_format::VIL_PIXEL_FORMAT_RGBA_UINT_16)
    .value("VIL_PIXEL_FORMAT_RGBA_INT_16", vil_pixel_format::VIL_PIXEL_FORMAT_RGBA_INT_16)
    .value("VIL_PIXEL_FORMAT_RGBA_BYTE", vil_pixel_format::VIL_PIXEL_FORMAT_RGBA_BYTE)
    .value("VIL_PIXEL_FORMAT_RGBA_SBYTE", vil_pixel_format::VIL_PIXEL_FORMAT_RGBA_SBYTE)
    .value("VIL_PIXEL_FORMAT_RGBA_FLOAT", vil_pixel_format::VIL_PIXEL_FORMAT_RGBA_FLOAT)
    .value("VIL_PIXEL_FORMAT_RGBA_DOUBLE", vil_pixel_format::VIL_PIXEL_FORMAT_RGBA_DOUBLE)
    //: std::complex<float> is a scalar for vil's purposes.
    .value("VIL_PIXEL_FORMAT_COMPLEX_FLOAT", vil_pixel_format::VIL_PIXEL_FORMAT_COMPLEX_FLOAT)
    //: std::complex<double> is a scalar for vil's purposes.
    .value("VIL_PIXEL_FORMAT_COMPLEX_DOUBLE", vil_pixel_format::VIL_PIXEL_FORMAT_COMPLEX_DOUBLE)
    .value("VIL_PIXEL_FORMAT_ENUM_END", vil_pixel_format::VIL_PIXEL_FORMAT_ENUM_END)
    // Add values here and be careful to keep values in vil_pixel_format.h/cxx in sync
    // Don't forget to increase the end value. Also add to vil_convert_cast in vil_convert.h
    .export_values();


  // Wrap vil_image_resource with a vil_smart_ptr handler
  // See: https://public.kitware.com/vxl/doc/development/books/core/book_7.html#SEC70
  // This is an exception. Normally, don't use holder types
  py::class_<vil_image_resource, PyImageResource /* <- trampoline */, vil_smart_ptr<vil_image_resource> /* <- holder type */ > (m, "image_resource")
    .def(py::init<>())
    .def("ni", &vil_image_resource::ni)
    .def("nj", &vil_image_resource::nj)
    .def("nplanes", &vil_image_resource::nplanes)
    .def("pixel_format", &vil_image_resource::pixel_format)
    .def("get_view_byte", [](vil_image_resource const& r) {vil_image_view<unsigned char> temp = r.get_view(); return temp; }, "Get a byte image view of all the data")
    .def("get_view_short", [](vil_image_resource const& r) {vil_image_view<unsigned short int> temp = r.get_view(); return temp; }, "Get a short image view of all the data")
    .def("get_view_float", [](vil_image_resource const& r) {vil_image_view<float> temp = r.get_view(); return temp; }, "Get a float image view of all the data")
    .def("get_view_int", [](vil_image_resource const& r) {vil_image_view<int> temp = r.get_view(); return temp; }, "Get an int image view of all the data")

    .def("get_view_byte", [](vil_image_resource const& r, unsigned i0, unsigned n_i, unsigned j0, unsigned n_j) {vil_image_view<unsigned char> temp = r.get_view(i0, n_i, j0, n_j); return temp; }, "Get a byte image view within a rectangular window")
    .def("get_view_short", [](vil_image_resource const& r, unsigned i0, unsigned n_i, unsigned j0, unsigned n_j) {vil_image_view<unsigned short int> temp = r.get_view(i0, n_i, j0, n_j); return temp; }, "Get a short image view within a rectangular window")
    .def("get_view_float", [](vil_image_resource const& r, unsigned i0, unsigned n_i, unsigned j0, unsigned n_j) {vil_image_view<float> temp = r.get_view(i0, n_i, j0, n_j); return temp; }, "Get a float image view within a rectangular window")
    .def("get_view_int", [](vil_image_resource const& r, unsigned i0, unsigned n_i, unsigned j0, unsigned n_j) {vil_image_view<int> temp = r.get_view(i0, n_i, j0, n_j); return temp; }, "Get an int image view within a rectangular window")

    .def("get_copy_view_byte", [](vil_image_resource const& r) {vil_image_view<unsigned char> temp = r.get_copy_view(); return temp; }, "Get a byte image view of a copy of all the data")
    .def("get_copy_view_short", [](vil_image_resource const& r) {vil_image_view<unsigned short int> temp = r.get_copy_view(); return temp; }, "Get a short image view of a copy of all the data")
    .def("get_copy_view_float", [](vil_image_resource const& r) {vil_image_view<float> temp = r.get_copy_view(); return temp; }, "Get a float image view of a copy of all the data")
    .def("get_copy_view_int", [](vil_image_resource const& r) {vil_image_view<int> temp = r.get_copy_view(); return temp; }, "Get an int image view of a copy of all the data")

    .def("get_copy_view_byte", [](vil_image_resource const& r, unsigned i0, unsigned n_i, unsigned j0, unsigned n_j) {vil_image_view<unsigned char> temp = r.get_copy_view(i0, n_i, j0, n_j); return temp; }, "Get a byte image view of a copy of this data within a rectangular window")
    .def("get_copy_view_short", [](vil_image_resource const& r, unsigned i0, unsigned n_i, unsigned j0, unsigned n_j) {vil_image_view<unsigned short int> temp = r.get_copy_view(i0, n_i, j0, n_j); return temp; }, "Get a short image view of a copy of this data within a rectangular window")
    .def("get_copy_view_float", [](vil_image_resource const& r, unsigned i0, unsigned n_i, unsigned j0, unsigned n_j) {vil_image_view<float> temp = r.get_copy_view(i0, n_i, j0, n_j); return temp; }, "Get a float image view of a copy of this data within a rectangular window")
    .def("get_copy_view_int", [](vil_image_resource const& r, unsigned i0, unsigned n_i, unsigned j0, unsigned n_j) {vil_image_view<int> temp = r.get_copy_view(i0, n_i, j0, n_j); return temp; }, "Get an int image view of a copy of this data within a rectangular window")

    .def("put_view", (bool (vil_image_resource::*)(const vil_image_view_base& im, unsigned i0, unsigned j0)) &vil_image_resource::put_view, "Put the data in this view back into the image source")
    .def("put_view", (bool (vil_image_resource::*)(const vil_image_view_base& im)) &vil_image_resource::put_view, "Put the data in this view back into the image source at the origin")
    .def("view_fits", &vil_image_resource::view_fits)
    .def("file_format", &vil_image_resource::file_format)
    .def("get_property", &vil_image_resource::get_property)
    .def("__len__", &resource_len)
    .def_property_readonly("shape", &image_resource_shape);


  // Need to provide a python wrapping for the base class so that
  // the derived classes work correctly.
  py::class_<vil_image_view_base, PyBaseView /* <- trampoline */> (m, "image_view_base")
    .def(py::init<>())
    .def("ni", &vil_image_view_base::ni)
    .def("nj", &vil_image_view_base::nj)
    .def("nplanes", &vil_image_view_base::nplanes)
    .def("size", &vil_image_view_base::size)
    .def("set_size", (void (vil_image_view_base::*)(unsigned, unsigned)) &vil_image_view_base::set_size, "Set size of current planes to width x height")
    .def("set_size", (void (vil_image_view_base::*)(unsigned, unsigned, unsigned)) &vil_image_view_base::set_size, "Resize to width x height x n_planes")
    .def("print", &vil_image_view_base::print)
    .def("is_a", &vil_image_view_base::is_a)
    .def("pixel_format", &vil_image_view_base::pixel_format)
    .def("is_class", &vil_image_view_base::is_class);

  // image view classes
  wrap_vil_image_view<bool>(m, "image_view_bool");
  wrap_vil_image_view<unsigned char>(m, "image_view_byte");
  wrap_vil_image_view<unsigned short int>(m, "image_view_uint16");
  wrap_vil_image_view<float>(m, "image_view_float");
  wrap_vil_image_view<int>(m, "image_view_int");
  wrap_vil_image_view<vil_rgb<unsigned char> >(m, "image_view_rgb_byte");

  m.def("new_image_resource_of_view", &vil_new_image_resource_of_view,
        "construct resource from view", py::arg("view"));

  m.def("crop_image_resource", (vil_image_resource_sptr (*)(const vil_image_resource_sptr&, unsigned, unsigned, unsigned, unsigned)) &vil_crop,
        py::arg("image_resource"), py::arg("i0"), py::arg("ni"), py::arg("j0"), py::arg("nj"));

  m.def("img_sum", &vil_image_sum_wrapper<unsigned char>, "", py::arg("image"), py::arg("p") = 0);
  m.def("img_sum", &vil_image_sum_wrapper<unsigned short int>, "", py::arg("image"), py::arg("p") = 0);
  m.def("img_sum", &vil_image_sum_wrapper<float>, "", py::arg("image"), py::arg("p") = 0);
  m.def("img_sum", &vil_image_sum_wrapper<int>, "", py::arg("image"), py::arg("p") = 0);

  m.def("_load_byte", &load_byte);
  m.def("_load_short", &load_short);
  m.def("_load_float", &load_float);
  m.def("_load_int", &load_int);

  m.def("load_image_resource", &vil_load_image_resource_wrapper);

  m.def("save_image_view", &vil_save_wrapper<unsigned char>);
  m.def("save_image_view", &vil_save_wrapper<unsigned short int>);
  m.def("save_image_view", &vil_save_wrapper<float>);
  m.def("save_image_view", &vil_save_wrapper<int>);

  m.def("_stretch_image_to_byte", &vil_stretch_image_to_byte_wrapper<unsigned char>,
        py::call_guard<py::gil_scoped_release>());
  m.def("_stretch_image_to_byte", &vil_stretch_image_to_byte_wrapper<unsigned short int>,
        py::call_guard<py::gil_scoped_release>());
  m.def("_stretch_image_to_byte", &vil_stretch_image_to_byte_wrapper<float>,
        py::call_guard<py::gil_scoped_release>());
  m.def("_stretch_image_to_byte", &vil_stretch_image_to_byte_wrapper<int>,
        py::call_guard<py::gil_scoped_release>());

  m.def("_stretch_image_to_short", &vil_stretch_image_to_short_wrapper<unsigned char>,
        py::call_guard<py::gil_scoped_release>());
  m.def("_stretch_image_to_short", &vil_stretch_image_to_short_wrapper<unsigned short int>,
        py::call_guard<py::gil_scoped_release>());
  m.def("_stretch_image_to_short", &vil_stretch_image_to_short_wrapper<float>,
        py::call_guard<py::gil_scoped_release>());
  m.def("_stretch_image_to_short", &vil_stretch_image_to_short_wrapper<int>,
        py::call_guard<py::gil_scoped_release>());

  m.def("_stretch_image_to_float", &vil_stretch_image_to_float_wrapper<unsigned char>,
        py::call_guard<py::gil_scoped_release>());
  m.def("_stretch_image_to_float", &vil_stretch_image_to_float_wrapper<unsigned short int>,
        py::call_guard<py::gil_scoped_release>());
  m.def("_stretch_image_to_float", &vil_stretch_image_to_float_wrapper<float>,
        py::call_guard<py::gil_scoped_release>());
  m.def("_stretch_image_to_float", &vil_stretch_image_to_float_wrapper<int>,
        py::call_guard<py::gil_scoped_release>());

  m.def("truncate_image_range", &vil_math_truncate_range<unsigned char>,
        py::call_guard<py::gil_scoped_release>());
  m.def("truncate_image_range", &vil_math_truncate_range<unsigned short int>,
        py::call_guard<py::gil_scoped_release>());
  m.def("truncate_image_range", &vil_math_truncate_range<float>,
        py::call_guard<py::gil_scoped_release>());
  m.def("truncate_image_range", &vil_math_truncate_range<int>,
        py::call_guard<py::gil_scoped_release>());

  m.def("image_range", &vil_image_range_wrapper<unsigned char>,
        py::call_guard<py::gil_scoped_release>());
  m.def("image_range", &vil_image_range_wrapper<unsigned short int>,
        py::call_guard<py::gil_scoped_release>());
  m.def("image_range", &vil_image_range_wrapper<float>,
        py::call_guard<py::gil_scoped_release>());
  m.def("image_range", &vil_image_range_wrapper<int>,
        py::call_guard<py::gil_scoped_release>());

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


PYBIND11_MODULE(_vil, m)
{
  m.doc() =  "Python bindings for the VIL computer vision libraries";

  pyvxl::vil::wrap_vil(m);
}
