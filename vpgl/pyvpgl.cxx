#include "pyvpgl.h"

#include <vpgl/vpgl_camera.h>
#include <vpgl/vpgl_proj_camera.h>
#include <vpgl/vpgl_affine_camera.h>
#include <vpgl/vpgl_perspective_camera.h>
#include <vpgl/vpgl_rational_camera.h>
#include <vpgl/vpgl_local_rational_camera.h>
#include <vpgl/vpgl_RSM_camera.h>
#include <vpgl/vpgl_lvcs.h>
#include <vpgl/vpgl_utm.h>
#include <vpgl/vpgl_fundamental_matrix.h>
#include <vpgl/vpgl_affine_fundamental_matrix.h>
#include <vpgl/vpgl_tri_focal_tensor.h>
#include <vpgl/vpgl_affine_tri_focal_tensor.h>
#include <vgl/vgl_box_2d.h>
#include <vgl/vgl_box_3d.h>
#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_ray_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_vector_2d.h>
#include <vgl/vgl_homg_point_2d.h>
#include <vgl/vgl_homg_point_3d.h>
#include <vgl/vgl_homg_line_2d.h>

#include <vpgl/file_formats/vpgl_geo_camera.h>
#include <vpgl/file_formats/vpgl_nitf_rational_camera.h>
#include <vpgl/file_formats/vpgl_nitf_RSM_camera_extractor.h>

#include <vil/vil_load.h>
#include <vil/vil_image_resource.h>
#include <vil/vil_image_resource_sptr.h>

#include "vul/vul_file.h"

// io classes for py::pickle
#include <vpgl/io/vpgl_io_proj_camera.h>
#include <vpgl/io/vpgl_io_affine_camera.h>
#include <vpgl/io/vpgl_io_perspective_camera.h>
#include <vpgl/io/vpgl_io_rational_camera.h>

#include "../pyvxl_util.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>

namespace py = pybind11;

namespace pyvxl { namespace vpgl {

/* This is a "trampoline" helper class for the virtual vpgl_camera<double>
 * class, which redirects virtual calls back to Python */
class PyCameraDouble : public vpgl_camera<double> {
public:
  using vpgl_camera<double>::vpgl_camera; /* Inherit constructors */

  std::string type_name() const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD(
      std::string,          /* Return type */
      vpgl_camera<double>,  /* Parent class */
      type_name,            /* Name of function */
                            /* No arguments, the trailing comma
                             * is needed for some compilers */
    );
  }

  void project(const double x, const double y, const double z, double& u, double& v) const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      void,                 /* Return type */
      vpgl_camera<double>,  /* Parent class */
      project,              /* Name of function */
      x, y, z, u, v         /* Arguments */
    );
  }

  PyCameraDouble* clone(void) const override {
    /* Generate wrapping code that enables native function overloading */
    PYBIND11_OVERLOAD_PURE(
      PyCameraDouble*,      /* Return type */
      vpgl_camera<double>,  /* Parent class */
      clone,                /* Name of function */
                            /* No arguments, the trailing comma
                             * is needed for some compilers */
    );
  }

};


template<class T>
vgl_point_2d<double> vpgl_project_point(T const& cam, vgl_point_3d<double> const& pt)
{
  double u,v;
  cam.project(pt.x(), pt.y(), pt.z(), u, v);
  return vgl_point_2d<double>(u,v);
}

template<class T>
vgl_vector_2d<double> vpgl_project_vector(T const& cam, vgl_vector_3d<double> const& vec)
{
  double u,v;
  cam.project(vec.x(), vec.y(), vec.z(), u, v);
  return vgl_vector_2d<double>(u,v);
}

template<class T>
vgl_homg_point_2d<double> vpgl_project_homg_point(T const& cam, vgl_homg_point_3d<double> const& x)
{
  return cam.project(x);
}

template<class T>
py::array vpgl_project_buffer(T const& cam, py::buffer b){

    py::buffer_info info = b.request();
    if(info.format != py::format_descriptor<double>::value){
        throw std::runtime_error("Incompatible scalar type");
    }

    if(info.ndim != 2){
        throw std::runtime_error("Expecting a 2-dimensional array");
    }

    if(info.shape[1] != 3){
        throw std::runtime_error("Expecting an Nx3 array");
    }

    double const* data = static_cast<double const*>(info.ptr);
    const size_t nextRow = info.strides[0] / sizeof(double);
    const size_t nextCol = info.strides[1] / sizeof(double);

    py::array output = py::array(py::buffer_info(
                (void*)nullptr, /* Numpy allocates */
                sizeof(double), /* Size of one item */
                py::format_descriptor<double>::value, /* Buffer format */
                2, /* Number of dimensions */
                std::vector<size_t>({static_cast<size_t>(info.shape[0]), 2}), /* Number of elements in each dimension */
                std::vector<size_t>({2*sizeof(double), sizeof(double)}) /* Strides for each dimension */
            ));
    py::buffer_info out_info = output.request();
    double* out_data = static_cast<double*>(out_info.ptr);
    const size_t output_nextRow = out_info.strides[0] / sizeof(double);
    const size_t output_nextCol = out_info.strides[1] / sizeof(double);

    for(py::ssize_t i = 0; i < info.shape[0]; ++i, data += nextRow, out_data += output_nextRow){
        const double x = *data;
        const double y = *(data + nextCol);
        const double z = *(data + 2 * nextCol);
        double u;
        double v;
        cam.project(x, y, z, u, v);
        *out_data = u;
        *(out_data + output_nextCol) = v;
    }

    return output;
}

template<class T>
std::tuple<double,double> vpgl_project_xyz(T const& cam, double x, double y, double z)
{
  double u,v;
  cam.project(x,y,z,u,v);
  return std::make_tuple(u,v);
}

vpgl_local_rational_camera<double> correct_local_rational_camera(vpgl_local_rational_camera<double> const& cam,
                                                                 double gt_offset_u, double gt_offset_v, bool verbose)
{
  if (verbose) {
    std::cout << "vxl.vpgl.correct_rational_camera, LOCAL rational camera, (off_u, off_v) = ("
              << gt_offset_u << ", " << gt_offset_v << ")" << std::endl;
  }

  // create copy
  vpgl_local_rational_camera<double> cam_out_local_rational(cam);

  // get u-v offset
  double offset_u, offset_v;
  cam_out_local_rational.image_offset(offset_u,offset_v);

  // add correction to offset
  offset_u += gt_offset_u;
  offset_v += gt_offset_v;

  // set new u-v offset
  cam_out_local_rational.set_image_offset(offset_u,offset_v);

  return cam_out_local_rational;
}

vpgl_rational_camera<double> correct_rational_camera(vpgl_rational_camera<double> const& cam_rational,
                                                     double gt_offset_u, double gt_offset_v, bool verbose)
{
  if (verbose) {
    std::cout << "vxl.vpgl.correct_rational_camera, rational camera, (off_u, off_v) = ("
              << gt_offset_u << ", " << gt_offset_v << ")" << std::endl;
  }

  // create copy
  vpgl_rational_camera<double> cam_out_rational(cam_rational);

  // get u-v offset
  double offset_u, offset_v;
  cam_out_rational.image_offset(offset_u, offset_v);

  // add correction to offset
  offset_u += gt_offset_u;
  offset_v += gt_offset_v;

  // set new u-v offset
  cam_out_rational.set_image_offset(offset_u, offset_v);

  return cam_out_rational;
}

vpgl_local_rational_camera<double>
__create_local_rational_camera(const vpgl_rational_camera<double>& rat_cam,
                             const vpgl_lvcs& lvcs,
                             unsigned min_x, unsigned min_y)
{
  // calculate local camera offset from image bounding box
  double global_u, global_v, local_u, local_v;
  rat_cam.image_offset(global_u, global_v);
  local_u = global_u - (double)min_x;  // the image was cropped by pixel
  local_v = global_v - (double)min_y;

  // create the local camera
  vpgl_local_rational_camera<double> local_camera(lvcs, rat_cam);
  local_camera.set_image_offset(local_u, local_v);
  return local_camera;
}

bool __project_box(const vpgl_rational_camera<double>& rat_cam, vpgl_lvcs &lvcs,
    const vgl_box_3d<double> &scene_bbox, double uncertainty,
    vgl_box_2d<double> &roi_box_2d)
{
  // project box
  double xoff, yoff, zoff;
  xoff = rat_cam.offset(vpgl_rational_camera<double>::X_INDX);
  yoff = rat_cam.offset(vpgl_rational_camera<double>::Y_INDX);
  zoff = rat_cam.offset(vpgl_rational_camera<double>::Z_INDX);

  // global to local (wgs84 to meter in order to apply uncertainty)
  double lx, ly, lz;
  lvcs.global_to_local(xoff, yoff, zoff, vpgl_lvcs::wgs84, lx, ly, lz, vpgl_lvcs::DEG, vpgl_lvcs::METERS);
  double center[3];
  center[0] = lx;  center[1] = ly;  center[2] = lz;

  // create a camera box with uncertainty
  vgl_box_3d<double> cam_box(center, 2*uncertainty, 2*uncertainty, 2*uncertainty, vgl_box_3d<double>::centre);
  std::vector<vgl_point_3d<double> > cam_corners = cam_box.vertices();

  // create the 3D box given input coordinates (in geo-coordinates)
  std::vector<vgl_point_3d<double> > box_corners = scene_bbox.vertices();

  // projection
  double lon, lat, gz;
  for (auto & cam_corner : cam_corners)
  {
    lvcs.local_to_global(cam_corner.x(), cam_corner.y(), cam_corner.z(), vpgl_lvcs::wgs84,
                          lon, lat, gz, vpgl_lvcs::DEG, vpgl_lvcs::METERS);
    vpgl_rational_camera<double>* new_cam = rat_cam.clone();
    new_cam->set_offset(vpgl_rational_camera<double>::X_INDX, lon);
    new_cam->set_offset(vpgl_rational_camera<double>::Y_INDX, lat);
    new_cam->set_offset(vpgl_rational_camera<double>::Z_INDX, gz);

    // project the box to image coords
    for (auto & box_corner : box_corners) {
      vgl_point_2d<double> p2d = new_cam->project(vgl_point_3d<double>(box_corner.x(), box_corner.y(), box_corner.z()));
      roi_box_2d.add(p2d);
    }
    delete new_cam;
  }

  return true;
}

double __clamp(double val, double min_val, double max_val)
{
  if (min_val > max_val) {
    throw std::runtime_error("invalid clip");
  }

  if (val < min_val)
    val = min_val;
  if (val > max_val)
    val = max_val;

  return val;
}

std::tuple<vpgl_local_rational_camera<double>, unsigned, unsigned, unsigned, unsigned>
crop_image_using_3d_box(
    unsigned img_ncols, unsigned img_nrows, vpgl_rational_camera<double> const& cam,
    double lower_left_lon, double lower_left_lat, double lower_left_elev,
    double upper_right_lon, double upper_right_lat, double upper_right_elev,
    double uncertainty, vpgl_lvcs lvcs)
{

  // generate a lvcs coordinates to transfer camera offset coordinates
  double ori_lon, ori_lat, ori_elev;
  lvcs.get_origin(ori_lat, ori_lon, ori_elev);
  if ( (ori_lat + ori_lon + ori_elev) * (ori_lat + ori_lon + ori_elev) < 1E-7) {
    lvcs = vpgl_lvcs(lower_left_lat, lower_left_lon, lower_left_elev, vpgl_lvcs::wgs84, vpgl_lvcs::DEG, vpgl_lvcs::METERS);
  }

  vgl_box_3d<double> scene_bbox(lower_left_lon, lower_left_lat, lower_left_elev,
                        upper_right_lon, upper_right_lat, upper_right_elev);

  vgl_box_2d<double> roi_box_2d;
  bool good = __project_box(cam, lvcs, scene_bbox, uncertainty, roi_box_2d);
  if(!good) {
    throw std::runtime_error("vxl.vpgl.crop_image_using_3d_box: failed to project 2d roi box");
  }
  std::cout << "vxl.vpgl.crop_image_using_3d_box: projected 2d roi box: " << roi_box_2d << " given uncertainty " << uncertainty << " meters." << std::endl;

  // convert to image coordinates
  auto x0 = (unsigned)__clamp(std::round(roi_box_2d.min_x()), 0, img_ncols-1);
  auto x1 = (unsigned)__clamp(std::round(roi_box_2d.max_x()), 0, img_ncols-1);
  auto y0 = (unsigned)__clamp(std::round(roi_box_2d.min_y()), 0, img_nrows-1);
  auto y1 = (unsigned)__clamp(std::round(roi_box_2d.max_y()), 0, img_nrows-1);

  // width & height
  auto nx = x1 - x0 + 1;
  auto ny = y1 - y0 + 1;

  if (nx <= 1 || ny <= 1)
  {
    throw std::runtime_error("vxl.vpgl.crop_image_using_3d_box: clipping box is out of image boundary, empty crop image returned");
  }

  // create the local camera
  vpgl_local_rational_camera<double> local_camera = __create_local_rational_camera(cam, lvcs, x0, y0);

  return std::make_tuple(local_camera, x0, y0, nx, ny);
}


void save_rational_camera(vpgl_camera<double> & cam, std::string camera_filename)
{
  auto *local_cam = dynamic_cast<vpgl_local_rational_camera<double>*>(&cam);

  if (local_cam) {
    if (!local_cam->save(camera_filename)) {
      std::ostringstream buffer;
      buffer << "Failed to save local rational camera " << camera_filename << std::endl;
      throw std::runtime_error(buffer.str());
    }
  }
  else {

    auto *rat_cam = dynamic_cast<vpgl_rational_camera<double>*>(&cam);

    if (!rat_cam) {
      throw std::runtime_error("error: could not convert camera input to a vpgl_rational_camera or local rational camera");
    }

    if (!rat_cam->save(camera_filename)) {
      std::ostringstream buffer;
      buffer << "Failed to save rational camera " << camera_filename << std::endl;
      throw std::runtime_error(buffer.str());
    }
  }
}

template<class CAM_T>
CAM_T load_camera(std::string const& camera_filename)
{
  // open file
  std::ifstream ifs(camera_filename.c_str());
  if (!ifs.is_open()) {
    std::ostringstream buffer;
    buffer << "Failed to open perspective camera file " << camera_filename << std::endl;
    throw std::runtime_error(buffer.str());
  }

  // create camera
  CAM_T cam;

  // load the file data into the camera
  std::string ext = vul_file_extension(camera_filename);
  if (ext == ".vsl") // binary form
  {
    vsl_b_ifstream bp_in(camera_filename.c_str());
    vsl_b_read(bp_in, cam);
    bp_in.close();
  }
  else {
   ifs >> cam;
  }
  return cam;
}

std::unique_ptr<vpgl_geo_camera> create_geocam(vnl_matrix<double> const& trans_matrix)
{
  if ((trans_matrix.rows() != 4) || (trans_matrix.cols() != 4)) {
    throw::std::runtime_error("trans_matrix should be of shape 4x4");
  }
  vpgl_lvcs_sptr lvcs(nullptr);  // no lvcs
  std::unique_ptr<vpgl_geo_camera> geocam(new vpgl_geo_camera(trans_matrix, lvcs));
  if ((trans_matrix(0,0) != 1.0) || (trans_matrix(1,1) != 1.0)) {
    geocam->set_scale_format(true);
  }
  return geocam;
}

std::unique_ptr<vpgl_geo_camera> create_geocam_with_lvcs(vnl_matrix<double> const& trans_matrix,
                                                         vpgl_lvcs const& lvcs)
{
  if ((trans_matrix.rows() != 4) || (trans_matrix.cols() != 4)) {
    throw::std::runtime_error("trans_matrix should be of shape 4x4");
  }
  vpgl_lvcs_sptr lvcs_copy(new vpgl_lvcs(lvcs));
  std::unique_ptr<vpgl_geo_camera> geocam(new vpgl_geo_camera(trans_matrix, lvcs_copy));
  if ((trans_matrix(0,0) != 1.0) || (trans_matrix(1,1) != 1.0)) {
    geocam->set_scale_format(true);
  }
  return geocam;
}

// vpgl_tri_focal_tensor "trampoline" helper class for virtual methods
template<class T>
class PyFundamentalMatrix : public vpgl_fundamental_matrix<T> {
 public:

  // Inherit constructors
  using vpgl_fundamental_matrix<T>::vpgl_fundamental_matrix;

  // set_matrix
  void set_matrix( const vnl_matrix_fixed<T,3,3>& F )
  {
    PYBIND11_OVERRIDE(void, vpgl_fundamental_matrix<T>,
                      set_matrix, F);
  }

};

// templated vpgl_fundamental_matrix
template<class T>
void wrap_vpgl_fundamental_matrix(py::module &m, const char* name)
{
  using VPGL_FM = vpgl_fundamental_matrix<T>;

  py::class_<VPGL_FM, PyFundamentalMatrix<T> /* <- trampoline */> (m, name)

    .def(py::init<>())
    .def(py::init<vnl_matrix_fixed<T,3,3> >())
    .def(py::init<vpgl_proj_camera<T>, vpgl_proj_camera<T> >(),
        py::arg("cr"), py::arg("cl"))

    // .def(py::init<vpgl_calibration_matrix<T>,
    //               vpgl_calibration_matrix<T>,
    //               vpgl_essential_matrix<T> >(),
    //     py::arg("kr"), py::arg("kl"), py::arg("em"))

    .def("__str__", streamToString<VPGL_FM>)

    .def_property_readonly("epipoles",
        [](VPGL_FM& self) {
          vgl_homg_point_2d<double> e12, e13;
          self.get_epipoles(e12, e13);
          return py::make_tuple(e12,e13);
        })

    .def("r_epipolar_line",
        overload_cast_<vgl_homg_point_2d<T> const&>
                      ()(&VPGL_FM::r_epipolar_line, py::const_),
        py::arg("pl"))
    .def("l_epipolar_line",
        overload_cast_<vgl_homg_point_2d<T> const&>
                      ()(&VPGL_FM::l_epipolar_line, py::const_),
        py::arg("pr"))

    .def("r_epipolar_line",
        overload_cast_<vgl_homg_line_2d<T> const&>
                      ()(&VPGL_FM::r_epipolar_line, py::const_),
        py::arg("epiline_l"))
    .def("l_epipolar_line",
        overload_cast_<vgl_homg_line_2d<T> const&>
                      ()(&VPGL_FM::l_epipolar_line, py::const_),
        py::arg("epiline_r"))

    .def("extract_left_camera",
        overload_cast_<vnl_vector_fixed<T,3> const&, T>
                      ()(&VPGL_FM::extract_left_camera, py::const_),
        py::arg("v"), py::arg("lambda"))
    // .def("extract_left_camera",
    //     overload_cast_<std::vector< vgl_point_3d<T> > const&,
    //                    std::vector< vgl_point_3d<T> > const&>
    //                   ()(&VPGL_FM::extract_left_camera, py::const_),
    //     py::arg("world_points"), py::arg("image_points"))

    .def("get_matrix", &VPGL_FM::get_matrix)
    .def("set_matrix",
        overload_cast_<vnl_matrix_fixed<T,3,3> const&>
                      ()(&VPGL_FM::set_matrix))
    .def("set_matrix",
        overload_cast_<vpgl_proj_camera<T> const&,
                       vpgl_proj_camera<T> const&>
                      ()(&VPGL_FM::set_matrix),
        py::arg("cr"), py::arg("cl"))

    .def_property("matrix",
        &VPGL_FM::get_matrix,
        overload_cast_<vnl_matrix_fixed<T,3,3> const&>
                      ()(&VPGL_FM::set_matrix))

    .def_property_readonly("svd", &VPGL_FM::svd)

    ;
}

// templated vpgl_affine_fundamental_matrix
template<class T>
void wrap_vpgl_affine_fundamental_matrix(py::module &m, const char* name)
{
  using VPGL_AFM = vpgl_affine_fundamental_matrix<T>;

  py::class_<VPGL_AFM, vpgl_fundamental_matrix<T> /* <- Parent */> (m, name)

    .def(py::init<>())
    .def(py::init<vnl_matrix_fixed<T,3,3> >())
    .def(py::init<vpgl_affine_camera<T>, vpgl_affine_camera<T> >(),
        py::arg("Ar"), py::arg("Al"))
    // vpgl_affine_fundamental_matrix( const vpgl_fundamental_matrix<T>& fm );

    .def("__str__", streamToString<VPGL_AFM>)

    .def("set_from_params", &VPGL_AFM::set_from_params,
        py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"), py::arg("e"))

    ;
}

// vpgl_tri_focal_tensor "trampoline" helper class for virtual methods
template<class T>
class PyTriFocalTensor : public vpgl_tri_focal_tensor<T> {
 public:

  // Inherit constructors
  using vpgl_tri_focal_tensor<T>::vpgl_tri_focal_tensor;

  // compute
  bool compute() override
  {
    PYBIND11_OVERRIDE(bool, vpgl_tri_focal_tensor<T>,
                      compute, );
  }

  // comparison
  // virtual bool operator==(vpgl_tri_focal_tensor<T> const& tensor) const

  // point transfer
  vgl_homg_point_2d<T> image1_transfer(
      vgl_homg_point_2d<T> const& point2,
      vgl_homg_point_2d<T> const& point3) const override
  {
    PYBIND11_OVERRIDE(vgl_homg_point_2d<T>, vpgl_tri_focal_tensor<T>,
                      image1_transfer, point2, point3);
  }

  vgl_homg_point_2d<T> image2_transfer(
      vgl_homg_point_2d<T> const& point1,
      vgl_homg_point_2d<T> const& point3) const override
  {
    PYBIND11_OVERRIDE(vgl_homg_point_2d<T>, vpgl_tri_focal_tensor<T>,
                      image2_transfer, point1, point3);
  }

  vgl_homg_point_2d<T> image3_transfer(
      vgl_homg_point_2d<T> const& point1,
      vgl_homg_point_2d<T> const& point2) const override
  {
    PYBIND11_OVERRIDE(vgl_homg_point_2d<T>, vpgl_tri_focal_tensor<T>,
                      image3_transfer, point1, point2);
  }

  // line transfer
  vgl_homg_line_2d<T> image1_transfer(
      vgl_homg_line_2d<T> const& line2,
      vgl_homg_line_2d<T> const& line3) const override
  {
    PYBIND11_OVERRIDE(vgl_homg_line_2d<T>, vpgl_tri_focal_tensor<T>,
                      image1_transfer, line2, line3);
  }

  vgl_homg_line_2d<T> image2_transfer(
      vgl_homg_line_2d<T> const& line1,
      vgl_homg_line_2d<T> const& line3) const override
  {
    PYBIND11_OVERRIDE(vgl_homg_line_2d<T>, vpgl_tri_focal_tensor<T>,
                      image2_transfer, line1, line3);
  }

  vgl_homg_line_2d<T> image3_transfer(
      vgl_homg_line_2d<T> const& line1,
      vgl_homg_line_2d<T> const& line2) const override
  {
    PYBIND11_OVERRIDE(vgl_homg_line_2d<T>, vpgl_tri_focal_tensor<T>,
                      image3_transfer, line1, line2);
  }

  // homographies induced by a line
  vgl_h_matrix_2d<T> hmatrix_13(vgl_homg_line_2d<T> const& line2) const override
  {
    PYBIND11_OVERRIDE(vgl_h_matrix_2d<T>, vpgl_tri_focal_tensor<T>,
                      hmatrix_13, line2);
  }

  vgl_h_matrix_2d<T> hmatrix_12(vgl_homg_line_2d<T> const& line3) const override
  {
    PYBIND11_OVERRIDE(vgl_h_matrix_2d<T>, vpgl_tri_focal_tensor<T>,
                      hmatrix_12, line3);
  }

  // epipoles
  void get_epipoles(vgl_homg_point_2d<T>& e12, vgl_homg_point_2d<T>& e13) override
  {
    PYBIND11_OVERRIDE(void, vpgl_tri_focal_tensor<T>,
                      get_epipoles, e12, e13);
  }

  vgl_homg_point_2d<T> epipole_12() override
  {
    PYBIND11_OVERRIDE(vgl_homg_point_2d<T>, vpgl_tri_focal_tensor<T>,
                      epipole_12, );
  }

  vgl_homg_point_2d<T> epipole_13() override
  {
    PYBIND11_OVERRIDE(vgl_homg_point_2d<T>, vpgl_tri_focal_tensor<T>,
                      epipole_13, );
  }

};


// templated vpgl_tri_focal_tensor
template<class T>
void wrap_vpgl_tri_focal_tensor(py::module &m, const char* name)
{
  using VPGL_TFT = vpgl_tri_focal_tensor<T>;

  py::class_<VPGL_TFT, PyTriFocalTensor<T> /* <- trampoline */> (m, name)

    .def(py::init<>())
    // vpgl_tri_focal_tensor(const vbl_array_3d<Type>&)
    // vpgl_tri_focal_tensor(const Type *tri_focal_tensor_array)
    .def(py::init<vpgl_proj_camera<T>, vpgl_proj_camera<T>, vpgl_proj_camera<T> >(),
        py::arg("camera1"), py::arg("camera2"), py::arg("camera3"))
    .def(py::init<vpgl_proj_camera<T>, vpgl_proj_camera<T> >(),
        py::arg("camera2"), py::arg("camera3"))
    .def(py::init<vnl_matrix_fixed<T,3,4>, vnl_matrix_fixed<T,3,4>, vnl_matrix_fixed<T,3,4> >(),
        py::arg("matrix1"), py::arg("matrix2"), py::arg("matrix3"))
    .def(py::init<vnl_matrix_fixed<T,3,4>, vnl_matrix_fixed<T,3,4> >(),
        py::arg("matrix2"), py::arg("matrix3"))
    .def("__str__", streamToString<VPGL_TFT>)
    .def(py::self == py::self)

    // index operators
    // T& operator()
    // T operator() const
    // void set(size_t i1, size_t i2, size_t i3, Type value)

    // setters
    .def("set",
        overload_cast_<vpgl_proj_camera<T> const&,
                       vpgl_proj_camera<T> const&,
                       vpgl_proj_camera<T> const&>
                      ()(&VPGL_TFT::set),
         py::arg("camera1"), py::arg("camera2"), py::arg("camera3"))
    .def("set",
        overload_cast_<vpgl_proj_camera<T> const&,
                       vpgl_proj_camera<T> const&>
                      ()(&VPGL_TFT::set),
         py::arg("camera2"), py::arg("camera3"))
    .def("set",
        overload_cast_<vnl_matrix_fixed<T,3,4> const&,
                       vnl_matrix_fixed<T,3,4> const&,
                       vnl_matrix_fixed<T,3,4> const&>
                      ()(&VPGL_TFT::set),
         py::arg("matrix1"), py::arg("matrix2"), py::arg("matrix3"))
    .def("set",
        overload_cast_<vnl_matrix_fixed<T,3,4> const&,
                       vnl_matrix_fixed<T,3,4> const&>
                      ()(&VPGL_TFT::set),
        py::arg("matrix2"), py::arg("matrix3"))

    // compute operations
    .def("compute", &VPGL_TFT::compute)
    .def("compute_epipoles", &VPGL_TFT::compute_epipoles)
    .def("compute_f_matrices", &VPGL_TFT::compute_f_matrices)
    .def("compute_proj_cameras", &VPGL_TFT::compute_proj_cameras)
    .def("compute_f_matrix_23", &VPGL_TFT::compute_f_matrix_23)

    // constraints
    .def("point_constraint_3x3", &VPGL_TFT::point_constraint_3x3,
        py::arg("point1"), py::arg("point2"), py::arg("point3"))
    .def("point_constraint", &VPGL_TFT::point_constraint,
        py::arg("point1"), py::arg("point2"), py::arg("point3"))
    .def("line_constraint_3", &VPGL_TFT::line_constraint_3,
        py::arg("line1"), py::arg("line2"), py::arg("line3"))

    // transfer
    .def("image1_point_transfer",
        overload_cast_<vgl_homg_point_2d<T> const&,
                       vgl_homg_point_2d<T> const&>
                      ()(&VPGL_TFT::image1_transfer, py::const_),
        py::arg("point2"), py::arg("point3"))
    .def("image2_point_transfer",
        overload_cast_<vgl_homg_point_2d<T> const&,
                       vgl_homg_point_2d<T> const&>
                      ()(&VPGL_TFT::image2_transfer, py::const_),
        py::arg("point1"), py::arg("point3"))
    .def("image3_point_transfer",
        overload_cast_<vgl_homg_point_2d<T> const&,
                       vgl_homg_point_2d<T> const&>
                      ()(&VPGL_TFT::image3_transfer, py::const_),
        py::arg("point1"), py::arg("point2"))

    .def("image1_line_transfer",
        overload_cast_<vgl_homg_line_2d<T> const&,
                       vgl_homg_line_2d<T> const&>
                      ()(&VPGL_TFT::image1_transfer, py::const_),
        py::arg("line2"), py::arg("line3"))
    .def("image2_line_transfer",
        overload_cast_<vgl_homg_line_2d<T> const&,
                       vgl_homg_line_2d<T> const&>
                      ()(&VPGL_TFT::image2_transfer, py::const_),
        py::arg("line1"), py::arg("line3"))
    .def("image3_line_transfer",
        overload_cast_<vgl_homg_line_2d<T> const&,
                       vgl_homg_line_2d<T> const&>
                      ()(&VPGL_TFT::image3_transfer, py::const_),
        py::arg("line1"), py::arg("line2"))

    // homographies induced by a line
    .def("hmatrix_13", &VPGL_TFT::hmatrix_13,
        py::arg("line2"))
    .def("hmatrix_12", &VPGL_TFT::hmatrix_12,
        py::arg("line3"))

    // epipoles
    .def_property_readonly("epipoles",
        [](VPGL_TFT& self) {
          vgl_homg_point_2d<double> e12, e13;
          self.get_epipoles(e12, e13);
          return py::make_tuple(e12,e13);
        })
    .def_property_readonly("epipole_12", &VPGL_TFT::epipole_12)
    .def_property_readonly("epipole_13", &VPGL_TFT::epipole_13)

    // fundamental matrices
    .def_property_readonly("fmatrix_12", &VPGL_TFT::fmatrix_12)
    .def_property_readonly("fmatrix_13", &VPGL_TFT::fmatrix_13)
    .def_property_readonly("fmatrix_23", &VPGL_TFT::fmatrix_23)

    // cameras
    .def_property_readonly("proj_camera_1", &VPGL_TFT::proj_camera_1)
    .def_property_readonly("proj_camera_2", &VPGL_TFT::proj_camera_2)
    .def_property_readonly("proj_camera_3", &VPGL_TFT::proj_camera_3)

    // constraint lines
    // void get_constraint_lines_image1(vgl_homg_point_2d<Type> const& p2,
    //                                  vgl_homg_point_2d<Type> const& p3,
    //                                  std::vector<vgl_homg_line_2d<Type> >& lines) const;

    // void get_constraint_lines_image2(vgl_homg_point_2d<Type> const& p1,
    //                                  vgl_homg_point_2d<Type> const& p3,
    //                                  std::vector<vgl_homg_line_2d<Type> >& lines) const;

    // void get_constraint_lines_image3(vgl_homg_point_2d<Type> const& p1,
    //                                  vgl_homg_point_2d<Type> const& p2,
    //                                  std::vector<vgl_homg_line_2d<Type> >& lines) const;

    // utility methods
    .def("postmultiply", &VPGL_TFT::postmultiply,
        py::arg("axis"), py::arg("matrix"))
    .def("premultiply", &VPGL_TFT::premultiply,
        py::arg("axis"), py::arg("matrix"))

    .def("postmultiply1", &VPGL_TFT::postmultiply1)
    .def("postmultiply2", &VPGL_TFT::postmultiply2)
    .def("postmultiply3", &VPGL_TFT::postmultiply3)

    .def("dot1", &VPGL_TFT::dot1)
    .def("dot2", &VPGL_TFT::dot2)
    .def("dot3", &VPGL_TFT::dot3)
    .def("dot1t", &VPGL_TFT::dot1t)
    .def("dot2t", &VPGL_TFT::dot2t)
    .def("dot3t", &VPGL_TFT::dot3t)

    ;
}


// templated vpgl_affine_tri_focal_tensor
template<class T>
void wrap_vpgl_affine_tri_focal_tensor(py::module &m, const char* name)
{
  using VPGL_ATFT = vpgl_affine_tri_focal_tensor<T>;

  py::class_<VPGL_ATFT, vpgl_tri_focal_tensor<T> /* <- Parent */> (m, name)

    .def(py::init<>())
    // vpgl_affine_tri_focal_tensor(const vbl_array_3d<T>&)
    // vpgl_affine_tri_focal_tensor(const vpgl_tri_focal_tensor<T>&)
    // vpgl_affine_tri_focal_tensor(const T *affine_tri_focal_tensor_array):
    .def(py::init<vpgl_affine_camera<T>, vpgl_affine_camera<T>, vpgl_affine_camera<T> >(),
        py::arg("camera1"), py::arg("camera2"), py::arg("camera3"))
    .def(py::init<vpgl_affine_camera<T>, vpgl_affine_camera<T> >(),
        py::arg("camera2"), py::arg("camera3"))
    .def(py::init<vnl_matrix_fixed<T,2,4>, vnl_matrix_fixed<T,2,4>, vnl_matrix_fixed<T,2,4> >(),
        py::arg("matrix1"), py::arg("matrix2"), py::arg("matrix3"))
    .def(py::init<vnl_matrix_fixed<T,2,4>, vnl_matrix_fixed<T,2,4> >(),
        py::arg("matrix2"), py::arg("matrix3"))

    // vpgl_affine_tri_focal_tensor(
    //     const vpgl_affine_camera<Type> &c1,
    //     const vpgl_affine_camera<Type> &c2,
    //     const vpgl_affine_camera<Type> &c3,
    //     std::vector<vgl_h_matrix_2d<Type>> img_pt_transforms)

    // vpgl_affine_tri_focal_tensor(
    //     const vpgl_affine_camera<Type>& c1,
    //     const vpgl_affine_camera<Type>& c2,
    //     const vpgl_affine_camera<Type>& c3,
    //     std::vector<std::pair<size_t, size_t> > const& image_dims_ni_nj)

    .def("__str__", streamToString<VPGL_ATFT>)

    // setters
    .def("set",
        overload_cast_<vpgl_affine_camera<T> const&,
                       vpgl_affine_camera<T> const&,
                       vpgl_affine_camera<T> const&>
                      ()(&VPGL_ATFT::set),
        py::arg("camera1"), py::arg("camera2"), py::arg("camera3"))
    .def("set",
        overload_cast_<vpgl_affine_camera<T> const&,
                       vpgl_affine_camera<T> const&>
                      ()(&VPGL_ATFT::set),
         py::arg("camera2"), py::arg("camera3"))
    .def("set",
        overload_cast_<vnl_matrix_fixed<T,2,4> const&,
                       vnl_matrix_fixed<T,2,4> const&,
                       vnl_matrix_fixed<T,2,4> const&>
                      ()(&VPGL_ATFT::set),
         py::arg("matrix1"), py::arg("matrix2"), py::arg("matrix3"))
    .def("set",
        overload_cast_<vnl_matrix_fixed<T,2,4> const&,
                       vnl_matrix_fixed<T,2,4> const&>
                      ()(&VPGL_ATFT::set),
         py::arg("matrix2"), py::arg("matrix3"))

    // set(const vpgl_affine_camera<Type> & c1,
    //     const vpgl_affine_camera<Type> & c2,
    //     const vpgl_affine_camera<Type> & c3,
    //     std::vector<vgl_h_matrix_2d<Type> > img_pt_transforms);

    // set(const vpgl_affine_camera<Type> & c1,
    //     const vpgl_affine_camera<Type> & c2,
    //     const vpgl_affine_camera<Type> & c3,
    //     std::vector<std::pair<size_t, size_t> > const & image_dims_ni_nj);

    // fundamental matrices
    .def_property_readonly("affine_fmatrix_12", &VPGL_ATFT::affine_fmatrix_12)
    .def_property_readonly("affine_fmatrix_13", &VPGL_ATFT::affine_fmatrix_13)
    .def_property_readonly("affine_fmatrix_23", &VPGL_ATFT::affine_fmatrix_23)

    // cameras
    .def_property_readonly("affine_camera_1", &VPGL_ATFT::affine_camera_1)
    .def_property_readonly("affine_camera_2", &VPGL_ATFT::affine_camera_2)
    .def_property_readonly("affine_camera_3", &VPGL_ATFT::affine_camera_3)

    ;
}

struct double3 {
  std::array<double, 3> value;
};

double3 lvcs_global_to_local_wrapper(
    vpgl_lvcs & lvcs, double lon, double lat, double el,
    vpgl_lvcs::cs_names const ics, vpgl_lvcs::AngUnits const iau,
    vpgl_lvcs::LenUnits const ilu
  )
{
    double3 result;
    lvcs.global_to_local(lon, lat, el, ics, result.value[0], result.value[1], result.value[2], iau, ilu);
    return result;
}

double3 lvcs_local_to_global_wrapper(
    vpgl_lvcs & lvcs, double x, double y, double z,
    vpgl_lvcs::cs_names const ocs, vpgl_lvcs::AngUnits const oau,
    vpgl_lvcs::LenUnits const olu
  )
{
    double3 result;
    lvcs.local_to_global(x, y, z, ocs, result.value[0], result.value[1], result.value[2], oau, olu);
    return result;
}

void wrap_vpgl(py::module &m)
{

  // This abstract class is the base for all the following cameras
  py::class_<vpgl_camera<double>, PyCameraDouble /* <- trampoline */> (m, "camera")
    .def(py::init<>())
    .def("type_name", &vpgl_camera<double>::type_name)
    .def("project", &vpgl_camera<double>::project)
    .def("is_a", &vpgl_camera<double>::is_a)
    .def("is_class", &vpgl_camera<double>::is_class)
    .def("clone", &vpgl_camera<double>::clone)
    .def("copy", &vpgl_camera<double>::clone)
    ;


  // =====PROJECTIVE CAMERA=====
  py::class_<vpgl_proj_camera<double>, vpgl_camera<double> /* <- Parent */> (m, "proj_camera")
    .def(py::init<vnl_matrix_fixed<double,3,4> >())
    .def("__str__", streamToString<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_homg_point<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_point<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_vector<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_xyz<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_buffer<vpgl_proj_camera<double> >)
    .def("get_matrix", &vpgl_proj_camera<double>::get_matrix, py::return_value_policy::copy)
    .def("set_matrix",
        overload_cast_<vnl_matrix_fixed<double,3,4> const&>
                      ()(&vpgl_proj_camera<double>::set_matrix))
    .def(py::self == py::self)
    .def(py::pickle(&vslPickleGetState<vpgl_proj_camera<double> >,
                    &vslPickleSetState<vpgl_proj_camera<double> >))
    .def("save", &vpgl_proj_camera<double>::save,
         "save camera to text file", py::arg("camera_filename"))
    ;

  m.def("load_proj_camera", &load_camera<vpgl_proj_camera<double> >,
        py::arg("camera_filename"));


  // =====AFFINE CAMERA=====
  py::class_<vpgl_affine_camera<double>, vpgl_proj_camera<double> /* <- Parent */> (m, "affine_camera")
    .def(py::init<vnl_matrix_fixed<double,3,4> >())
    .def(py::init<vgl_vector_3d<double>, vgl_vector_3d<double>, vgl_point_3d<double>,
         double, double, double, double>(),
         py::arg("ray"), py::arg("up"), py::arg("stare_pt"),
         py::arg("u0"), py::arg("v0"), py::arg("su"), py::arg("sv"))
    .def("backproject_ray",
      [](vpgl_affine_camera<double> &cam, double u, double v){
        vgl_homg_point_2d<double> image_point(u, v);
        vgl_ray_3d<double> ray = cam.backproject_ray(image_point);
        return ray;
      }
      )
    .def("ray_dir", &vpgl_affine_camera<double>::ray_dir)
    .def_property("viewing_distance",
                  &vpgl_affine_camera<double>::viewing_distance,  // getter
                  &vpgl_affine_camera<double>::set_viewing_distance)  // setter
    /* .def("set_viewing_distance", &vpgl_affine_camera<double>::set_viewing_distance); */
    .def(py::pickle(&vslPickleGetState<vpgl_affine_camera<double> >,
                    &vslPickleSetState<vpgl_affine_camera<double> >))
    ;

  m.def("load_affine_camera", &load_camera<vpgl_affine_camera<double> >,
        py::arg("camera_filename"));


  // =====PERSPECTIVE CAMERA=====
  py::class_<vpgl_calibration_matrix<double> >(m, "calibration_matrix")
    .def(py::init<vnl_matrix_fixed<double,3,3> >())
    .def(py::init<double, vgl_point_2d<double> >())
    .def("get_matrix",&vpgl_calibration_matrix<double>::get_matrix)
    .def(py::self == py::self)
    ;

  py::class_<vpgl_perspective_camera<double>, vpgl_proj_camera<double> /* <- Parent */ > (m, "perspective_camera")
    .def(py::init<vpgl_calibration_matrix<double>, vgl_rotation_3d<double>, vgl_vector_3d<double> >())
    .def(py::init<vpgl_calibration_matrix<double>, vgl_point_3d<double>, vgl_rotation_3d<double> >())
    .def("__str__", streamToString<vpgl_perspective_camera<double> >)
    .def_property("camera_center",
                  &vpgl_perspective_camera<double>::get_camera_center,
                  &vpgl_perspective_camera<double>::set_camera_center)
    .def_property("calibration",
                  &vpgl_perspective_camera<double>::get_calibration,
                  &vpgl_perspective_camera<double>::set_calibration)
    .def_property("rotation",
                  &vpgl_perspective_camera<double>::get_rotation,
                  &vpgl_perspective_camera<double>::set_rotation)
    .def_property("translation",
                  &vpgl_perspective_camera<double>::get_translation,
                  &vpgl_perspective_camera<double>::set_translation)
    .def("principal_axis", &vpgl_perspective_camera<double>::principal_axis,
         "compute the principal axis (i.e. the vector perpendicular to the image plane pointing towards the front of the camera")
    .def("is_behind_camera", &vpgl_perspective_camera<double>::is_behind_camera,
         "Determine whether the given homogeneous world point lies in front of the principal plane",
         py::arg("world_point"))
    .def("backproject",
      [](vpgl_perspective_camera<double> &cam, double u, double v){
        vgl_line_3d_2_points<double> line2pts = cam.backproject(u,v);
        return line2pts;
      }
      )
    .def("backproject_ray",
      [](vpgl_perspective_camera<double> &cam, double u, double v){
        vgl_homg_point_2d<double> image_point(u, v);
        vgl_ray_3d<double> ray = cam.backproject_ray(image_point);
        return ray;
      }
      )
    .def(py::pickle(&vslPickleGetState<vpgl_perspective_camera<double> >,
                    &vslPickleSetState<vpgl_perspective_camera<double> >))
    ;

  m.def("load_perspective_camera", &load_camera<vpgl_perspective_camera<double> >,
        py::arg("camera_filename"));


  // =====SCALE OFFSET=====
  py::class_<vpgl_scale_offset<double> >(m, "scale_offset")
    .def(py::init<>())
    .def(py::init<double, double>())
    .def("set_scale", &vpgl_scale_offset<double>::set_scale)
    .def("set_offset", &vpgl_scale_offset<double>::set_offset)
    .def("scale", &vpgl_scale_offset<double>::scale)
    .def("offset", &vpgl_scale_offset<double>::offset)
    .def("normalize", &vpgl_scale_offset<double>::normalize)
    .def("un_normalize", &vpgl_scale_offset<double>::un_normalize)
    .def("__str__", [](const vpgl_scale_offset<double>& scale_offset){
        std::ostringstream buffer;
        buffer << "< scale_offset<double> scale = " << scale_offset.scale();
        buffer << ", offset = " << scale_offset.offset() << " >";
        return buffer.str();
    });


  // =====RATIONAL CAMERA=====
  py::class_<vpgl_rational_camera<double>, vpgl_camera<double> /* <- Parent */ > rational_camera(m, "rational_camera");

  // enumerations, attached to parent class
  py::enum_<vpgl_rational_order> rational_order(m, "rational_order");
  for (auto item : vpgl_rational_order_func::initializer_list) {
    rational_order.value(vpgl_rational_order_func::to_string(item).c_str(), item);
  }

  // enumerations, attached to this class
  py::enum_<vpgl_rational_camera<double>::coor_index>(rational_camera, "coor_index")
    .value("X_INDX", vpgl_rational_camera<double>::X_INDX)
    .value("Y_INDX", vpgl_rational_camera<double>::Y_INDX)
    .value("Z_INDX", vpgl_rational_camera<double>::Z_INDX)
    .value("U_INDX", vpgl_rational_camera<double>::U_INDX)
    .value("V_INDX", vpgl_rational_camera<double>::V_INDX);

  py::enum_<vpgl_rational_camera<double>::poly_index>(rational_camera, "poly_index")
    .value("NEU_U", vpgl_rational_camera<double>::NEU_U)
    .value("DEN_U", vpgl_rational_camera<double>::DEN_U)
    .value("NEU_V", vpgl_rational_camera<double>::NEU_V)
    .value("DEN_V", vpgl_rational_camera<double>::DEN_V);

  // function definitions
  rational_camera

    // overloaded constructors
    .def(py::init<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>,
                  double, double, double, double, double, double, double, double, double, double,
                  vpgl_rational_order>(),
        py::arg("neu_u"), py::arg("den_u"), py::arg("neu_v"), py::arg("den_v"),
        py::arg("x_scale"), py::arg("x_off"), py::arg("y_scale"), py::arg("y_off"),
        py::arg("z_scale"), py::arg("z_off"), py::arg("u_scale"), py::arg("u_off"),
        py::arg("v_scale"), py::arg("v_off"),
        py::arg("input_rational_order") = vpgl_rational_order::VXL)

    .def(py::init<vnl_matrix_fixed<double,4,20>, std::vector<vpgl_scale_offset<double> >,
                  vpgl_rational_order>(),
        py::arg("rational_coeffs"), py::arg("scale_offsets"),
        py::arg("input_rational_order") = vpgl_rational_order::VXL)

    // general functions
    .def("__str__", streamToString<vpgl_rational_camera<double> >)
    .def(py::self == py::self)
    .def(py::pickle(&vslPickleGetState<vpgl_rational_camera<double> >,
                    &vslPickleSetState<vpgl_rational_camera<double> >))

    // save to file
    .def("save", &vpgl_rational_camera<double>::save,
        py::arg("cam_path"),
        py::arg("rational_order") = vpgl_rational_order::RPC00B)

    // point projection
    .def("project", vpgl_project_point<vpgl_rational_camera<double> >)
    .def("project", vpgl_project_buffer<vpgl_rational_camera<double> >)
    .def("project", vpgl_project_xyz<vpgl_rational_camera<double> >)

    // getter/setter
    .def("coefficient_matrix", &vpgl_rational_camera<double>::coefficient_matrix,
        py::arg("rational_order") = vpgl_rational_order::VXL)
    .def("scale_offsets", &vpgl_rational_camera<double>::scale_offsets)
    .def("offset", &vpgl_rational_camera<double>::offset)
    .def("scale", &vpgl_rational_camera<double>::scale)
    .def_property("image_offset",
        [](vpgl_rational_camera<double>& self) {
          double u,v; self.image_offset(u,v);
          return py::make_tuple(u,v);
        },
        [](vpgl_rational_camera<double>& self, const std::array<double,2>& uv) {
          self.set_image_offset(uv[0],uv[1]);
        })
    .def_property("image_scale",
        [](vpgl_rational_camera<double>& self) {
          double u,v; self.image_scale(u,v);
          return py::make_tuple(u,v);
        },
        [](vpgl_rational_camera<double>& self, const std::array<double,2>& uv) {
          self.set_image_scale(uv[0],uv[1]);
        })
     ;

  // Functions related to rational cameras
  m.def("_read_rational_camera",
        [](std::string const& fname){return read_rational_camera<double>(fname);},
        py::return_value_policy::take_ownership);
  m.def("_read_rational_camera_from_txt",
        [](std::string const& fname){return read_rational_camera_from_txt<double>(fname);},
        py::return_value_policy::take_ownership);
  m.def("load_rational_camera_from_str",
        [](std::string const& input){
          std::istringstream iss(input);
          return read_rational_camera<double>(iss);
        },
        py::return_value_policy::take_ownership);
  m.def("_correct_rational_camera", &correct_rational_camera);
  m.def("save_rational_camera", &save_rational_camera,
        py::arg("cam"), py::arg("camera_filename"));


  py::class_<vpgl_local_rational_camera<double>, vpgl_rational_camera<double> /* <- Parent */ > (m, "local_rational_camera")
    .def(py::init<vpgl_lvcs const&, vpgl_rational_camera<double> const&>())
    .def(py::init<double, double, double, vpgl_rational_camera<double> const&>())
    .def("set_lvcs",
         (void (vpgl_local_rational_camera<double>::*)(vpgl_lvcs const&))
             &vpgl_local_rational_camera<double>::set_lvcs,
         py::arg("lvcs"))
    .def("set_lvcs",
         (void (vpgl_local_rational_camera<double>::*)(double const&, double const&, double const&))
             &vpgl_local_rational_camera<double>::set_lvcs,
         py::arg("longitude"),py::arg("latitude"),py::arg("elevation"))
    .def("lvcs", &vpgl_local_rational_camera<double>::lvcs)
    .def("clone_as_rational_camera",
         [](py::object &self) { return self.cast<vpgl_rational_camera<double> >(); },
         "Return a copy of the camera as a vpgl_rational_camera")
    ;

  m.def("_read_local_rational_camera",
        [](std::string const& fname){return read_local_rational_camera<double>(fname);},
        py::return_value_policy::take_ownership);
  m.def("_correct_local_rational_camera", &correct_local_rational_camera);


  // =====NITF RATIONAL CAMERA=====
  py::class_<vpgl_nitf_rational_camera, vpgl_rational_camera<double> /* <- Parent */ > nitf_rational_camera(m, "nitf_rational_camera");

  // enumerations, attached to this class
  py::enum_<vpgl_nitf_rational_camera::geopt_coord>(nitf_rational_camera, "geopt_coord")
    .value("LON", vpgl_nitf_rational_camera::LON)
    .value("LAT", vpgl_nitf_rational_camera::LAT);

  py::enum_<vpgl_nitf_rational_camera::igeolo_order>(nitf_rational_camera, "igeolo_order")
    .value("UL", vpgl_nitf_rational_camera::UL)
    .value("UR", vpgl_nitf_rational_camera::UR)
    .value("LR", vpgl_nitf_rational_camera::LR)
    .value("LL", vpgl_nitf_rational_camera::LL);

  // function definitions
  nitf_rational_camera

    // constructors
    .def(py::init<>())
    .def(py::init<std::string const&, bool>(),
         py::arg("file"), py::arg("verbose")=false)

    // properties
    .def_property_readonly("rational_extension_type", &vpgl_nitf_rational_camera::rational_extension_type)
    .def_property_readonly("image_id", &vpgl_nitf_rational_camera::image_id)

    // read from file
    .def("read",
         overload_cast_<std::string const&, bool>()(&vpgl_nitf_rational_camera::read),
         py::arg("file"), py::arg("verbose")=false)
    ;

  // load rational camera from nitf file
  m.def("load_rational_camera_nitf",
        [] (std::string const& file, bool verbose) {
          auto camera = vpgl_nitf_rational_camera();
          if (!camera.read(file, verbose)) {
            throw std::runtime_error("Unable to load RPCs from NITF");
          }
          return camera;
        },
        py::arg("file"), py::arg("verbose")=false);


  // ===== REPLACEMENT SENSOR MODEL =====

  // enumeration
  py::enum_<vpgl_ground_domain_id>(m, "ground_domain_id")
    .value("G", vpgl_ground_domain_id::G)
    .value("H", vpgl_ground_domain_id::H)
    .value("R", vpgl_ground_domain_id::R)
    ;

  // ground domain
  py::class_<vpgl_ground_domain<double> > (m, "ground_domain")
    .def(py::init<>())
    .def(py::init<vpgl_ground_domain_id>(), py::arg("id"))
    .def(py::init<std::string>(), py::arg("id"))

    .def("__str__", streamToString<vpgl_ground_domain<double> >)
    .def("__repr__", streamToString<vpgl_ground_domain<double> >)
    .def("as_dict", struct_to_dict<vpgl_ground_domain<double> >)

    .def_readwrite("id", &vpgl_ground_domain<double>::id_)
    .def_readwrite("translation", &vpgl_ground_domain<double>::translation_)
    .def_readwrite("rotation", &vpgl_ground_domain<double>::rotation_)

    .def("reset", &vpgl_ground_domain<double>::reset)

    .def("world_to_ground",
        [] (const vpgl_ground_domain<double> & self,
            const double lon, const double lat, const double elev)
        {
          double x, y, z;
          self.world_to_ground(lon, lat, elev, x, y, z);
          return std::make_tuple(x, y, z);
        },
        py::arg("lon"), py::arg("lat"), py::arg("elev"))

    .def("world_to_ground",
        overload_cast_<vnl_vector_fixed<double, 3> const&>
                      ()(&vpgl_ground_domain<double>::world_to_ground, py::const_),
        py::arg("world_point"))

    .def("world_to_ground",
        overload_cast_<vgl_point_3d<double> const&>
                      ()(&vpgl_ground_domain<double>::world_to_ground, py::const_),
        py::arg("world_point"))

    ;

  // polycam
  py::class_<vpgl_polycam<double> > polycam(m, "polycam");

  // enumerations, attached to this class
  py::enum_<vpgl_polycam<double>::coor_index>(polycam, "coor_index")
    .value("X_INDX", vpgl_polycam<double>::X_INDX)
    .value("Y_INDX", vpgl_polycam<double>::Y_INDX)
    .value("Z_INDX", vpgl_polycam<double>::Z_INDX)
    .value("U_INDX", vpgl_polycam<double>::U_INDX)
    .value("V_INDX", vpgl_polycam<double>::V_INDX);

  py::enum_<vpgl_polycam<double>::poly_index>(polycam, "poly_index")
    .value("NEU_U", vpgl_polycam<double>::NEU_U)
    .value("DEN_U", vpgl_polycam<double>::DEN_U)
    .value("NEU_V", vpgl_polycam<double>::NEU_V)
    .value("DEN_V", vpgl_polycam<double>::DEN_V);

  py::enum_<vpgl_polycam<double>::poly_comp_index>(polycam, "poly_comp_index")
    .value("NEU_U", vpgl_polycam<double>::P_NEU_U)
    .value("DEN_U", vpgl_polycam<double>::P_DEN_U)
    .value("NEU_V", vpgl_polycam<double>::P_NEU_V)
    .value("DEN_V", vpgl_polycam<double>::P_DEN_V);

  // function definitions
  polycam

    // constructors
    .def(py::init<>())
    .def(py::init<size_t, size_t>(),
        py::arg("ridx"), py::arg("cidx"))
    .def(py::init<size_t, size_t,
                  const std::vector<std::vector<int>> &,
                  const std::vector<std::vector<double>> &,
                  const std::vector<vpgl_scale_offset<double>> &>(),
        py::arg("ridx"), py::arg("cidx"),
        py::arg("powers"), py::arg("coefficients"), py::arg("scale_offsets"))

    // accessors
    .def_property_readonly("ridx", &vpgl_polycam<double>::ridx)
    .def_property_readonly("cidx", &vpgl_polycam<double>::cidx)
    ;

  // RSM Camera
  py::class_<vpgl_RSM_camera<double>, vpgl_camera<double> /* <- Parent */ > RSM_camera(m, "RSM_camera");

  // function definitions
  RSM_camera

    // overloaded constructors
    // default only since construction is complex
    // and is done by the vpgl_nitf_RSM_extractor factory class
    .def(py::init<>())
    .def(py::init<const vpgl_polycam<double> &>(),
        py::arg("polycam"))
    .def(py::init(
        [](const vpgl_polycam<double> & pcam,
           const vpgl_ground_domain<double> & gd)
        {
          auto rsm_camera = vpgl_RSM_camera<double>(pcam);
          rsm_camera.set_ground_domain(gd);
          return rsm_camera;
        }),
        py::arg("polycam"), py::arg("ground_domain"))
    // point projection
    .def("project", vpgl_project_point<vpgl_RSM_camera<double> >)
    .def("project", vpgl_project_buffer<vpgl_RSM_camera<double> >)
    .def("project", vpgl_project_xyz<vpgl_RSM_camera<double> >)
    // adjustable parameters image row (v)  and column (u) offsets
    // apply to the entire image even if it has multiple polynomial sections
    .def("set_adjustable_parameters", &vpgl_RSM_camera<double>::set_adjustable_parameters, py::arg("adj_u"),py::arg("adj_v"))
    .def("adjustable_parameters", &vpgl_RSM_camera<double>::adjustable_parameters)
    .def("n_regions", &vpgl_RSM_camera<double>::n_regions, " number of polynomial regions")
    .def_property("ground_domain",
                  &vpgl_RSM_camera<double>::ground_domain,
                  &vpgl_RSM_camera<double>::set_ground_domain)
    ;

  // =====LOCAL VERTICAL COORDINATE SYSTEM (LVCS)=====
  py::class_<vpgl_lvcs> lvcs(m, "lvcs");


  // enumerations, attached to LVCS class
  py::enum_<vpgl_lvcs::LenUnits>(lvcs, "LenUnits")
    .value("FEET", vpgl_lvcs::FEET)
    .value("METERS", vpgl_lvcs::METERS);

  py::enum_<vpgl_lvcs::AngUnits>(lvcs, "AngUnits")
    .value("RADIANS", vpgl_lvcs::RADIANS)
    .value("DEG", vpgl_lvcs::DEG);

  py::enum_<vpgl_lvcs::cs_names>(lvcs, "cs_names")
    .value("wgs84", vpgl_lvcs::wgs84)
    .value("nad27n", vpgl_lvcs::nad27n)
    .value("wgs72", vpgl_lvcs::wgs72)
    .value("utm", vpgl_lvcs::utm)
    .value("NumNames", vpgl_lvcs::NumNames);

  // function definitions
  lvcs

    // overloaded constructors
    .def(py::init<double,double,double,vpgl_lvcs::cs_names,double,double,vpgl_lvcs::AngUnits,vpgl_lvcs::LenUnits,double,double,double>(),
        py::arg("orig_lat")=0,py::arg("orig_lon")=0,py::arg("orig_elev")=0,
        py::arg("cs_name")=vpgl_lvcs::wgs84,
        py::arg("lat_scale")=0,py::arg("lon_scale")=0,
        py::arg("ang_unit")=vpgl_lvcs::DEG,py::arg("elev_unit")=vpgl_lvcs::METERS,
        py::arg("lox")=0,py::arg("loy")=0,py::arg("theta")=0)

    .def(py::init<double,double,double,vpgl_lvcs::cs_names,vpgl_lvcs::AngUnits,vpgl_lvcs::LenUnits>(),
        py::arg("orig_lat"),py::arg("orig_lon"),py::arg("orig_elev"),
        py::arg("cs_name")=vpgl_lvcs::wgs84,py::arg("ang_unit")=vpgl_lvcs::DEG,py::arg("elev_unit")=vpgl_lvcs::METERS)

    .def(py::init<double,double,double,double,double,vpgl_lvcs::cs_names,vpgl_lvcs::AngUnits,vpgl_lvcs::LenUnits>(),
        py::arg("lat_low"),py::arg("lon_low"),py::arg("lat_high"),py::arg("lon_high"),py::arg("elev"),
        py::arg("cs_name")=vpgl_lvcs::wgs84,py::arg("ang_unit")=vpgl_lvcs::DEG,py::arg("elev_unit")=vpgl_lvcs::METERS)

    // python print
    .def("__str__", streamToString<vpgl_lvcs>)
    .def("__repr__", streamToString<vpgl_lvcs>)

    // getters
    .def("get_origin",     [](vpgl_lvcs &L) {double lon,lat,e; L.get_origin(lat,lon,e); return std::make_tuple(lon,lat,e); })
    .def("get_scale",      [](vpgl_lvcs &L) {double lon,lat; L.get_scale(lat,lon); return std::make_tuple(lon,lat); })
    .def("get_transform",  [](vpgl_lvcs &L) {double lox,loy,th; L.get_transform(lox,loy,th); return std::make_tuple(lox,loy,th); })
    .def("get_utm_origin", [](vpgl_lvcs &L) {double x,y,e; int z; L.get_utm_origin(x,y,e,z); return std::make_tuple(x,y,e,z); })
    .def("get_utm",        [](vpgl_lvcs &L) {int z; bool sf; L.get_utm(z,sf); return std::make_tuple(z,sf); })
    .def("get_cs_name",    &vpgl_lvcs::get_cs_name)
    .def("get_len_unit",   &vpgl_lvcs::local_length_unit)
    .def("get_ang_unit",   &vpgl_lvcs::geo_angle_unit)

    // setters
    .def("set_transform", &vpgl_lvcs::set_transform,
         py::arg("lox"), py::arg("loy"), py::arg("theta"))
    .def("set_origin", &vpgl_lvcs::set_origin,
         py::arg("lon"), py::arg("lat"), py::arg("elev"))
    .def("set_utm", &vpgl_lvcs::set_utm,
         py::arg("zone"), py::arg("south_flag"))

    // read/write to string
    .def("reads",
        [](vpgl_lvcs &L, std::string const &str)
        {
          std::istringstream iss(str.c_str());
          if (iss) {
            L.read(iss);
            return true;
          } else
            return false;
        }
      )

    .def("writes",
        [](vpgl_lvcs &L)
        {
          std::ostringstream oss;
          L.write(oss);
          return oss.str();
        }
      )


    // read/write to file
    .def("read",
        [](vpgl_lvcs &L, std::string const &filename)
        {
          std::ifstream ifs(filename.c_str());
          if (ifs) {
            L.read(ifs);
            ifs.close();
            return true;
          } else
            return false;
        }
      )

    .def("write",
        [](vpgl_lvcs &L, std::string const &filename)
        {
          std::ofstream ofs(filename.c_str());
          if (ofs) {
            L.write(ofs); ofs.close();
            return true;
          } else
            return false;
        }
      )

    // local->global coordinate transform
    .def("_local_to_global",
        [](vpgl_lvcs &L, double const lx, double const ly, double const lz,
           vpgl_lvcs::cs_names const ocs, vpgl_lvcs::AngUnits const oau,
           vpgl_lvcs::LenUnits const olu)
          {
            double glon, glat, gz;
            L.local_to_global(lx,ly,lz,ocs,glon,glat,gz,oau,olu);
            return std::make_tuple(glon,glat,gz);
          },
        py::arg("local_x"),py::arg("local_y"),py::arg("local_z"),
        py::arg("output_cs_name"),py::arg("output_ang_unit")=vpgl_lvcs::DEG,
        py::arg("output_len_unit")=vpgl_lvcs::METERS
     )
     .def("_local_to_global", py::vectorize(lvcs_local_to_global_wrapper),
        py::arg("local_x"),py::arg("local_y"),py::arg("local_z"),
        py::arg("output_cs_name"),py::arg("output_ang_unit")=vpgl_lvcs::DEG,
        py::arg("output_len_unit")=vpgl_lvcs::METERS
     )

    // global->local coordinate transform
    .def("_global_to_local",
        [](vpgl_lvcs &L, const double glon, const double glat, const double gz,
           vpgl_lvcs::cs_names const ics, vpgl_lvcs::AngUnits const iau,
           vpgl_lvcs::LenUnits const ilu)
          {
            double lx, ly, lz;
            L.global_to_local(glon,glat,gz,ics,lx,ly,lz,iau,ilu);
            return std::make_tuple(lx,ly,lz);
          },
        py::arg("global_longitude"),py::arg("global_latitude"),py::arg("global_elevation"),
        py::arg("input_cs_name"),py::arg("input_ang_unit")=vpgl_lvcs::DEG,
        py::arg("input_len_unit")=vpgl_lvcs::METERS
     )
     .def("_global_to_local", py::vectorize(lvcs_global_to_local_wrapper),
          py::arg("global_longitude"),py::arg("global_latitude"),py::arg("global_elevation"),
          py::arg("input_cs_name"),py::arg("input_ang_unit")=vpgl_lvcs::DEG,
          py::arg("input_len_unit")=vpgl_lvcs::METERS
     )
    ;

  m.def("create_lvcs",
       [](double lat, double lon, double elev, std::string name)
       {return new vpgl_lvcs(lat, lon, elev,
                             vpgl_lvcs::str_to_enum(name.c_str()), vpgl_lvcs::DEG, vpgl_lvcs::METERS);},
       py::return_value_policy::take_ownership);


  // =====LAT/LON to UTM CONVERTER=====
  py::class_<vpgl_utm>(m, "utm")
    .def(py::init<>())
    .def("lonlat2utm",
        [] (vpgl_utm &U, double lon, double lat)
          { double x,y; int zone; U.transform(lat,lon,x,y,zone); return std::make_tuple(x,y,zone); },
        py::arg("longitude"),py::arg("latitude"))
    .def("utm2lonlat",
        [] (vpgl_utm &U, double x, double y, int zone, bool is_south)
          { double lat,lon; U.transform(zone,x,y,lat,lon,is_south); return std::make_tuple(lon,lat); },
        py::arg("easting"),py::arg("northing"),py::arg("zone"),py::arg("is_south")=false)
    ;


  // =====GEO-CAMERA=====
  // Geo- Camera definitions
  py::class_<vpgl_geo_camera, vpgl_camera<double> /* <- Parent */ > (m, "geo_camera")
    // Default methods
    .def(py::init<>())
    .def(py::init(&create_geocam))
    .def(py::init(&create_geocam_with_lvcs))
    .def("__str__", streamToString<vpgl_geo_camera >)
    // Convert pixel coords (u,v) to a lon/lat pair
    .def("img_to_global",
      [](vpgl_geo_camera &G, double const u, double const v)
      {
        double lon, lat;
        G.img_to_global(u, v, lon, lat);
        return std::make_tuple(lon, lat);
      },
      py::arg("u"), py::arg("v")
    )
    .def("get_geotransform",
      [](vpgl_geo_camera &G)
      {
        /* return GeoTransform suitable for GDAL GeoTiff */
        /* Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2) */
        /* Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5) */
        /* http://www.gdal.org/gdal_datamodel.html */
        /* Note that the pixel/line coordinates in the above are from (0.0,0.0) at the top */
        /* left corner of the top left pixel to (width_in_pixels,height_in_pixels) at the */
        /* bottom right corner of the bottom right pixel. The pixel/line location of the */
        /* center of the top left pixel would therefore be (0.5,0.5). */
        vnl_matrix<double> trans_matrix = G.trans_matrix();
        return std::make_tuple((double)trans_matrix[0][3],   // GT(0)
                               (double)trans_matrix[0][0],   // GT(1)
                               (double)trans_matrix[1][0],   // GT(2)
                               (double)trans_matrix[1][3],   // GT(3)
                               (double)trans_matrix[0][1],   // GT(4)
                               (double)trans_matrix[1][1]);  // GT(5)
      })
    .def("global_to_img",
         [](vpgl_geo_camera &G, double const lon, double const lat, const double elev)
         {
           double u,v;
           G.global_to_img(lon, lat, elev, u, v);
           return std::make_tuple(u,v);
         },
         py::arg("longitude"), py::arg("latitude"), py::arg("elevation")
        );

  m.def("load_geo_camera_from_geotiff", &load_geo_camera_from_geotiff,
        "Load a vpgl_geo_camera from a geotiff file",
        py::arg("file"), py::arg("lvcs")=nullptr);

  m.def("load_geo_camera_from_resource", &load_geo_camera_from_resource,
        "Load a vpgl_geo_camera from an image resource",
        py::arg("resource"), py::arg("lvcs")=nullptr);

  m.def("load_geo_camera_from_geotransform", &load_geo_camera_from_geotransform,
        "Load a vpgl_geo_camera from a GDAL geotransform",
        py::arg("geotransform"), py::arg("utm_zone")=-1,
        py::arg("northing")=0, py::arg("lvcs")=nullptr);

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


  // =====MISC=====
  // image cropping extents from 3D box and rational camera
  // Note this must be listed after the lvcs wrapper to enable the default lvcs
  m.def("crop_image_using_3d_box", &crop_image_using_3d_box,
        "Determine image cropping extents given 3D box & rational camera",
        py::call_guard<py::gil_scoped_release>(),
        py::arg("img_ncols"), py::arg("img_nrows"), py::arg("cam"),
        py::arg("lower_left_lon"), py::arg("lower_left_lat"), py::arg("lower_left_elev"),
        py::arg("upper_right_lon"), py::arg("upper_right_lat"), py::arg("upper_right_elev"),
        py::arg("uncertainty"), py::arg("lvcs") = vpgl_lvcs());


  // =====TEMPLATED=====
  wrap_vpgl_fundamental_matrix<double>(m, "fundamental_matrix");
  wrap_vpgl_affine_fundamental_matrix<double>(m, "affine_fundamental_matrix");

  wrap_vpgl_tri_focal_tensor<double>(m, "tri_focal_tensor");
  wrap_vpgl_affine_tri_focal_tensor<double>(m, "affine_tri_focal_tensor");


  //===== Replacement Sensor Model related interfaces
  py::class_<rsm_metadata> (m, "rsm_metadata")
    .def(py::init<>())
    .def("as_dict", struct_to_dict<rsm_metadata>)
    .def_readonly("any_valid", &rsm_metadata::any_valid)
    .def_readonly("platform_name", &rsm_metadata::platform_name_)
    .def_readonly("platform_name_valid", &rsm_metadata::platform_name_valid)
    .def_readonly("image_iid1", &rsm_metadata::image_iid1_)
    .def_readonly("image_iid1_valid", &rsm_metadata::image_iid1_valid)
    .def_readonly("image_iid2", &rsm_metadata::image_iid2_)
    .def_readonly("image_iid2_valid", &rsm_metadata::image_iid2_valid)
    .def_readonly("sid", &rsm_metadata::sid_)
    .def_readonly("sid_valid", &rsm_metadata::sid_valid)
    .def_readonly("acquisition_time", &rsm_metadata::acquisition_time_)
    .def_readonly("acquisition_time_valid", &rsm_metadata::acquisition_time_valid)
    .def_readonly("image_type", &rsm_metadata::image_type_)
    .def_readonly("image_type_valid", &rsm_metadata::image_type_valid)
    .def_readonly("ground_domain_valid", &rsm_metadata::ground_domain_valid)
    .def_readonly("ground_domain", &rsm_metadata::ground_domain_)
    .def_readonly("local_transform_valid", &rsm_metadata::local_transform_valid)
    .def_readonly("translation", &rsm_metadata::translation_)
    .def_readonly("rotation", &rsm_metadata::rotation_)
    .def_readonly("igeolo", &rsm_metadata::igeolo_)
    .def_readonly("igeolo_valid", &rsm_metadata::igeolo_valid)
    .def_readonly("polytope_valid", &rsm_metadata::polytope_valid)
    .def_readonly("polytope", &rsm_metadata::polytope_)
    .def_readonly("bounding_box_valid", &rsm_metadata::bounding_box_valid)
    .def_readonly("bounding_box", &rsm_metadata::bounding_box_)
    .def_readonly("corners_valid", &rsm_metadata::xy_corners_valid)
    .def_readonly("upper_left", &rsm_metadata::upper_left_)
    .def_readonly("upper_right", &rsm_metadata::upper_right_)
    .def_readonly("lower_left", &rsm_metadata::lower_left_)
    .def_readonly("lower_right", &rsm_metadata::lower_right_)
    .def_readonly("footprint", &rsm_metadata::footprint_)
    .def_readonly("image_corners_valid", &rsm_metadata::image_corners_valid)
    .def_readonly("min_image_corner", &rsm_metadata::min_image_corner_)
    .def_readonly("max_image_corner", &rsm_metadata::max_image_corner_)
    .def_readonly("sun_azimuth_radians", &rsm_metadata::sun_azimuth_radians_)
    .def_readonly("sun_azimuth_valid", &rsm_metadata::sun_azimuth_valid)
    .def_readonly("sun_elevation_radians", &rsm_metadata::sun_elevation_radians_)
    .def_readonly("sun_elevation_valid", &rsm_metadata::sun_elevation_valid)
    ;


  py::class_<RSM_ECA_adjustable_parameter_metadata> (m, "RSM_ECA_adjustable_parameter_metadata")
    .def(py::init<>())
    .def("print", &RSM_ECA_adjustable_parameter_metadata::print, py::arg("os"))
    .def("as_dict", struct_to_dict<RSM_ECA_adjustable_parameter_metadata>)

    .def_readonly("defined", &RSM_ECA_adjustable_parameter_metadata::defined_)
    .def_readonly("num_adj_params", &RSM_ECA_adjustable_parameter_metadata::num_adj_params_)
    .def_readonly("num_original_adj_params", &RSM_ECA_adjustable_parameter_metadata::num_original_adj_params_)
    .def_readonly("num_independent_subgroups", &RSM_ECA_adjustable_parameter_metadata::num_independent_subgroups_)
    .def_readonly("unmodeled_error", &RSM_ECA_adjustable_parameter_metadata::unmodeled_error_)

    .def_readonly("translation", &RSM_ECA_adjustable_parameter_metadata::translation_)
    .def_readonly("rotation", &RSM_ECA_adjustable_parameter_metadata::rotation_)

    .def_readonly("covar_index", &RSM_ECA_adjustable_parameter_metadata::covar_index_)

    .def_readonly("independent_subgroup_covariance", &RSM_ECA_adjustable_parameter_metadata::independent_subgroup_covariance_)

    .def_readonly("correlation_flags", &RSM_ECA_adjustable_parameter_metadata::correlation_flags_)
    .def_readonly("correlation_segments", &RSM_ECA_adjustable_parameter_metadata::correlation_segments_)
    .def_readonly("phi", &RSM_ECA_adjustable_parameter_metadata::phi_)

    .def_readonly("unmodeled_row_variance", &RSM_ECA_adjustable_parameter_metadata::unmodeled_row_variance_)
    .def_readonly("unmodeled_col_variance", &RSM_ECA_adjustable_parameter_metadata::unmodeled_col_variance_)
    .def_readonly("unmodeled_row_col_variance", &RSM_ECA_adjustable_parameter_metadata::unmodeled_row_col_variance_)
    .def_readonly("unmodeled_row_correlation", &RSM_ECA_adjustable_parameter_metadata::unmodeled_row_correlation_)
    .def_readonly("unmodeled_col_correlation", &RSM_ECA_adjustable_parameter_metadata::unmodeled_col_correlation_)

    ;

  py::class_<RSM_ECB_adjustable_parameter_metadata> (m, "RSM_ECB_adjustable_parameter_metadata")
    .def(py::init<>())
    .def("print", &RSM_ECB_adjustable_parameter_metadata::print, py::arg("os"))
    .def("as_dict", struct_to_dict<RSM_ECB_adjustable_parameter_metadata>)

    .def_readonly("defined", &RSM_ECB_adjustable_parameter_metadata::defined_)
    .def_readonly("unmodeled_error", &RSM_ECB_adjustable_parameter_metadata::unmodeled_error_)
    .def_readonly("num_active_adj_params", &RSM_ECB_adjustable_parameter_metadata::num_active_adj_params_)
    .def_readonly("num_original_adj_params", &RSM_ECB_adjustable_parameter_metadata::num_original_adj_params_)
    .def_readonly("num_independent_subgroups", &RSM_ECB_adjustable_parameter_metadata::num_independent_subgroups_)
    .def_readonly("image_adjustable_params", &RSM_ECB_adjustable_parameter_metadata::image_adjustable_params_)
    .def_readonly("rect_local_coordinate_system", &RSM_ECB_adjustable_parameter_metadata::rect_local_coordinate_system_)

    .def_readonly("n_image_adjustable_params", &RSM_ECB_adjustable_parameter_metadata::n_image_adjustable_params_)
    .def_readonly("n_image_row_adjustable_params", &RSM_ECB_adjustable_parameter_metadata::n_image_row_adjustable_params_)
    .def_readonly("n_image_col_adjustable_params", &RSM_ECB_adjustable_parameter_metadata::n_image_col_adjustable_params_)

    .def_readonly("image_row_xyz_powers", &RSM_ECB_adjustable_parameter_metadata::image_row_xyz_powers_)
    .def_readonly("image_col_xyz_powers", &RSM_ECB_adjustable_parameter_metadata::image_col_xyz_powers_)
    .def_readonly("ground_adjustable_params", &RSM_ECB_adjustable_parameter_metadata::ground_adjustable_params_)
    .def_readonly("n_ground_adjustable_params", &RSM_ECB_adjustable_parameter_metadata::n_ground_adjustable_params_)
    .def_readonly("ground_adjust_param_idx", &RSM_ECB_adjustable_parameter_metadata::ground_adjust_param_idx_)

    .def_readonly("basis_option", &RSM_ECB_adjustable_parameter_metadata::basis_option_)
    .def_readonly("n_basis_adjustable_params", &RSM_ECB_adjustable_parameter_metadata::n_basis_adjustable_params_)
    .def_readonly("A_matrix", &RSM_ECB_adjustable_parameter_metadata::A_matrix_)

    .def_readonly("independent_covar", &RSM_ECB_adjustable_parameter_metadata::independent_covar_)
    .def_readonly("correlation_domain_flags", &RSM_ECB_adjustable_parameter_metadata::correlation_domain_flags_)
    .def_readonly("corr_analytic_functions", &RSM_ECB_adjustable_parameter_metadata::corr_analytic_functions_)
    .def_readonly("corr_piecewise_functions", &RSM_ECB_adjustable_parameter_metadata::corr_piecewise_functions_)

    .def_readonly("mapping_matrix", &RSM_ECB_adjustable_parameter_metadata::mapping_matrix_)

    .def_readonly("xyz_norm", &RSM_ECB_adjustable_parameter_metadata::xyz_norm_)

    .def_readonly("rect_translation", &RSM_ECB_adjustable_parameter_metadata::rect_translation_)
    .def_readonly("rect_rotation", &RSM_ECB_adjustable_parameter_metadata::rect_rotation_)

    .def_readonly("unmodeled_row_variance", &RSM_ECB_adjustable_parameter_metadata::unmodeled_row_variance_)
    .def_readonly("unmodeled_col_variance", &RSM_ECB_adjustable_parameter_metadata::unmodeled_col_variance_)
    .def_readonly("unmodeled_row_col_variance", &RSM_ECB_adjustable_parameter_metadata::unmodeled_row_col_variance_)

    .def_readonly("unmodeled_analytic", &RSM_ECB_adjustable_parameter_metadata::unmodeled_analytic_)
    .def_readonly("unmodeled_row_analytic_function", &RSM_ECB_adjustable_parameter_metadata::unmodeled_row_analytic_function_)
    .def_readonly("unmodeled_col_analytic_function", &RSM_ECB_adjustable_parameter_metadata::unmodeled_col_analytic_function_)
    .def_readonly("unmodeled_row_piecewise_function", &RSM_ECB_adjustable_parameter_metadata::unmodeled_row_piecewise_function_)
    .def_readonly("unmodeled_col_piecewise_function", &RSM_ECB_adjustable_parameter_metadata::unmodeled_col_piecewise_function_)
    ;

  py::class_<vpgl_nitf_RSM_camera_extractor> (m, "vpgl_nitf_RSM_camera_extractor")
    .def(py::init<>())
    .def(py::init<std::string const&>(),"construct from file", py::arg("nitf_image_path"))
    .def("image_id", &vpgl_nitf_RSM_camera_extractor::image_id, "get name from metadata", py::arg("image_subheader_index"))

    .def("nitf_header_contains_RSM_tres", &vpgl_nitf_RSM_camera_extractor::nitf_header_contains_RSM_tres, "number of headers containing RSM data, a return of zero indicates no RSM present")
    .def("process", &vpgl_nitf_RSM_camera_extractor::process, "extract header TREs", py::arg("verbose"))
    .def("first_index_with_RSM", &vpgl_nitf_RSM_camera_extractor::first_index_with_RSM, "with multiple image subheaders, the first containing RSM data")

    .def("RSM_camera", &vpgl_nitf_RSM_camera_extractor::RSM_camera, "the RSM camera at a given index", py::arg("image_subheader_index")=0)
    .def("meta", &vpgl_nitf_RSM_camera_extractor::meta, "general metadata including sun angle", py::arg("image_subheader_index")=0)
    .def("RSM_ECA_adjustable_parameter_data", &vpgl_nitf_RSM_camera_extractor::RSM_ECA_adjustable_parameter_data, py::arg("image_subheader_index") = 0)
    .def("RSM_ECB_adjustable_parameter_data", &vpgl_nitf_RSM_camera_extractor::RSM_ECB_adjustable_parameter_data, py::arg("image_subheader_index") = 0)
    .def("save_tre_values", &vpgl_nitf_RSM_camera_extractor::save_tre_values, py::arg("file"))
    ;

}


  }}

PYBIND11_MODULE(_vpgl, m)
{
  m.doc() =  "Python bindings for the VPGL computer vision libraries";

  pyvxl::vpgl::wrap_vpgl(m);
  // required for return types of lvcs _global_to_local and _local_to_global
  PYBIND11_NUMPY_DTYPE(pyvxl::vpgl::double3, value);
}
