#include "pyvpgl.h"

#include <vpgl/vpgl_camera.h>
#include <vpgl/vpgl_proj_camera.h>
#include <vpgl/vpgl_affine_camera.h>
#include <vpgl/vpgl_perspective_camera.h>
#include <vpgl/vpgl_rational_camera.h>
#include <vpgl/vpgl_local_rational_camera.h>
#include <vpgl/vpgl_lvcs.h>
#include <vpgl/vpgl_utm.h>
#include <vgl/vgl_box_2d.h>
#include <vgl/vgl_box_3d.h>
#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_vector_2d.h>
#include <vgl/vgl_homg_point_2d.h>
#include <vgl/vgl_homg_point_3d.h>

#include <vpgl/file_formats/vpgl_geo_camera.h>
#include <vpgl/file_formats/vpgl_nitf_rational_camera.h>

#include <vil/vil_load.h>
#include <vil/vil_image_resource.h>
#include <vil/vil_image_resource_sptr.h>

#include "../pyvxl_util.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
namespace py = pybind11;

namespace pyvxl {

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
};

/* /1* This is a "trampoline" helper class for the virtual methods in */
/*  * vpgl_rational_camera<double>, which redirects virtual calls back to Python *1/ */
/* class PyRationalCameraDouble : public vpgl_rational_camera<double> { */
/* public: */
/*   using vpgl_rational_camera<double>::vpgl_rational_camera; /1* Inherit constructors *1/ */

/*   std::string type_name() const override { */
/*     /1* Generate wrapping code that enables native function overloading *1/ */
/*     PYBIND11_OVERLOAD( */
/*       std::string,                   /1* Return type *1/ */
/*       vpgl_rational_camera<double>,  /1* Parent class *1/ */
/*       type_name,                     /1* Name of function *1/ */
/*                                      /1* No arguments, the trailing comma */
/*                                       * is needed for some compilers *1/ */
/*     ); */
/*   } */

/*   void project(const double x, const double y, const double z, double& u, double& v) const override { */
/*     /1* Generate wrapping code that enables native function overloading *1/ */
/*     PYBIND11_OVERLOAD_PURE( */
/*       void,                          /1* Return type *1/ */
/*       vpgl_rational_camera<double>,  /1* Parent class *1/ */
/*       project,                       /1* Name of function *1/ */
/*       x, y, z, u, v                  /1* Arguments *1/ */
/*     ); */
/*   } */

/*   void print(std::ostream& s=std::cout, vpgl_rational_order output_order=vpgl_rational_order::VXL) const override { */
/*     /1* Generate wrapping code that enables native function overloading *1/ */
/*     PYBIND11_OVERLOAD_PURE( */
/*       void,                          /1* Return type *1/ */
/*       vpgl_rational_camera<double>,  /1* Parent class *1/ */
/*       print,                         /1* Name of function *1/ */
/*       s, output_order                /1* Arguments *1/ */
/*     ); */
/*   } */

/*   bool save(std::string cam_path, vpgl_rational_order output_order=vpgl_rational_order::RPC00B) const override { */
/*     /1* Generate wrapping code that enables native function overloading *1/ */
/*     PYBIND11_OVERLOAD_PURE( */
/*       bool,                          /1* Return type *1/ */
/*       vpgl_rational_camera<double>,  /1* Parent class *1/ */
/*       save,                          /1* Name of function *1/ */
/*       cam_path, output_order         /1* Arguments *1/ */
/*     ); */
/*   } */

/*   void write_pvl(std::ostream& s, vpgl_rational_order output_order) const override { */
/*     /1* Generate wrapping code that enables native function overloading *1/ */
/*     PYBIND11_OVERLOAD_PURE( */
/*       void,                          /1* Return type *1/ */
/*       vpgl_rational_camera<double>,  /1* Parent class *1/ */
/*       write_pvl,                     /1* Name of function *1/ */
/*       s, output_order                /1* Arguments *1/ */
/*     ); */
/*   } */

/* }; */



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

    for(size_t i = 0; i < info.shape[0]; ++i, data += nextRow, out_data += output_nextRow){
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

void wrap_vpgl(py::module &m)
{

  // This abstract class is the base for all the following cameras
  py::class_<vpgl_camera<double>, PyCameraDouble /* <- trampoline */> (m, "camera")
    .def(py::init<>())
    .def("type_name", &vpgl_camera<double>::type_name)
    .def("project", &vpgl_camera<double>::project)
    .def("is_a", &vpgl_camera<double>::is_a)
    .def("is_class", &vpgl_camera<double>::is_class);

  py::class_<vpgl_proj_camera<double>, vpgl_camera<double> /* <- Parent */> (m, "proj_camera")
    .def(py::init<vnl_matrix_fixed<double,3,4> >())
    .def("__str__", streamToString<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_homg_point<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_point<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_vector<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_xyz<vpgl_proj_camera<double> >)
    .def("get_matrix", &vpgl_proj_camera<double>::get_matrix, py::return_value_policy::copy);

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
                  &vpgl_affine_camera<double>::set_viewing_distance);  // setter
    /* .def("set_viewing_distance", &vpgl_affine_camera<double>::set_viewing_distance); */

  py::class_<vpgl_calibration_matrix<double> >(m, "calibration_matrix")
    .def(py::init<vnl_matrix_fixed<double,3,3> >())
    .def(py::init<double, vgl_point_2d<double> >());

  py::class_<vpgl_perspective_camera<double>, vpgl_proj_camera<double> /* <- Parent */ > (m, "perspective_camera")
    .def(py::init<vpgl_calibration_matrix<double>, vgl_rotation_3d<double>, vgl_vector_3d<double> >())
    .def("__str__", streamToString<vpgl_perspective_camera<double> >);

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
  py::enum_<vpgl_rational_order>(m, "rational_order")
    .value("VXL", vpgl_rational_order::VXL)
    .value("RPC00B", vpgl_rational_order::RPC00B)
    ;

  // enumerations, attached to this class
  py::enum_<vpgl_rational_camera<double>::coor_index>(rational_camera, "coor_index")
    .value("X_INDX", vpgl_rational_camera<double>::X_INDX)
    .value("Y_INDX", vpgl_rational_camera<double>::Y_INDX)
    .value("Z_INDX", vpgl_rational_camera<double>::Z_INDX)
    .value("U_INDX", vpgl_rational_camera<double>::U_INDX)
    .value("V_INDX", vpgl_rational_camera<double>::V_INDX);

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

    // copy constructor
    .def("copy", &vpgl_rational_camera<double>::clone)

    // python print
    .def("__str__", streamToString<vpgl_rational_camera<double> >)

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
    ;

  m.def("_read_local_rational_camera",
        [](std::string const& fname){return read_local_rational_camera<double>(fname);},
        py::return_value_policy::take_ownership);
  m.def("_correct_local_rational_camera", &correct_local_rational_camera);


  py::class_<vpgl_nitf_rational_camera, vpgl_rational_camera<double> /* <- Parent */ > (m, "nitf_rational_camera");
    /* .def(py::init<>()) */
    /* .def(py::init<std::string const&, bool>()) */
    /* .def(py::init<vil_nitf2_image*, bool>()); */

  m.def("_read_rational_camera_nitf",
        [](std::string const& fname){return read_rational_camera_from_txt<double>(fname);},
        py::return_value_policy::take_ownership);



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
    .def("get_cs_name",    &vpgl_lvcs::get_cs_name)
    .def("get_len_unit",   &vpgl_lvcs::local_length_unit)
    .def("get_ang_unit",   &vpgl_lvcs::geo_angle_unit)

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
    .def("local_to_global",
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

    // global->local coordinate transform
    .def("global_to_local",
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
        py::arg("input_len_unit")=vpgl_lvcs::METERS)

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
  m.def("crop_image_using_3d_box", crop_image_using_3d_box,
        "Determine image cropping extents given 3D box & rational camera",
        py::arg("img_ncols"), py::arg("img_nrows"), py::arg("cam"),
        py::arg("lower_left_lon"), py::arg("lower_left_lat"), py::arg("lower_left_elev"),
        py::arg("upper_right_lon"), py::arg("upper_right_lat"), py::arg("upper_right_elev"),
        py::arg("uncertainty"), py::arg("lvcs") = vpgl_lvcs());

}
}

PYBIND11_MODULE(_vpgl, m)
{
  m.doc() =  "Python bindings for the VIL computer vision libraries";

  pyvxl::wrap_vpgl(m);
}
