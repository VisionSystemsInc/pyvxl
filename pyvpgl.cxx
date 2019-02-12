#include "pyvpgl.h"
#include <vpgl/vpgl_proj_camera.h>
#include <vpgl/vpgl_affine_camera.h>
#include <vpgl/vpgl_perspective_camera.h>
#include <vpgl/vpgl_rational_camera.h>
#include <vpgl/vpgl_local_rational_camera.h>
#include <vpgl/vpgl_lvcs.h>
#include <vpgl/vpgl_utm.h>
#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_homg_point_2d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_vector_2d.h>
#include <vgl/vgl_homg_point_2d.h>
#include <vgl/vgl_homg_point_3d.h>
#include "pyvxl_util.h"

#include <vpgl/file_formats/vpgl_geo_camera.h>
#include <vil/vil_load.h>
#include <vil/vil_image_resource.h>
#include <vil/vil_image_resource_sptr.h>

 // For rational to affine //
#include <vnl/vnl_random.h>
#include <vpgl/vpgl_local_rational_camera.h>
#include <vpgl/algo/vpgl_camera_compute.h>
#include <vpgl/algo/vpgl_affine_rectification.h>
#include <vgl/vgl_ray_3d.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <sstream>
#include <vector>

namespace py = pybind11;

namespace pyvxl {

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

void wrap_vpgl(py::module &m)
{
  py::class_<vpgl_proj_camera<double> >(m, "proj_camera")
    .def(py::init<vnl_matrix_fixed<double,3,4> >())
    .def("__str__", stream2str<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_homg_point<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_point<vpgl_proj_camera<double> >)
    .def("project", vpgl_project_vector<vpgl_proj_camera<double> >)
    .def("get_matrix", &vpgl_proj_camera<double>::get_matrix, py::return_value_policy::copy);

  py::class_<vpgl_affine_camera<double>, vpgl_proj_camera<double> >(m, "affine_camera")
    .def(py::init<vnl_matrix_fixed<double,3,4> >())
    .def("__str__", stream2str<vpgl_proj_camera<double> >)
    .def("project",
      [](vpgl_affine_camera<double> &cam, double x, double y, double z)
      {
        double u, v;
        cam.project(x, y, z, u, v);
        return std::make_tuple(u, v);
      }
    )
    .def("backproject_ray",
      [](vpgl_affine_camera<double> &cam, double u, double v){
        vgl_homg_point_2d<double> image_point(u, v);
        vgl_ray_3d<double> ray = cam.backproject_ray(image_point);
        return ray;
      }
    );

  py::class_<vpgl_calibration_matrix<double> >(m, "calibration_matrix")
    .def(py::init<vnl_matrix_fixed<double,3,3> >())
    .def(py::init<double, vgl_point_2d<double> >());

  py::class_<vpgl_perspective_camera<double>, vpgl_proj_camera<double> >(m, "perspective_camera")
    .def(py::init<vpgl_calibration_matrix<double>, vgl_rotation_3d<double>, vgl_vector_3d<double> >())
    .def("__str__", stream2str<vpgl_perspective_camera<double> >);

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

  py::class_<vpgl_rational_camera<double> > rational_camera(m, "rational_camera");
  py::enum_<vpgl_rational_camera<double>::coor_index>(rational_camera, "coor_index")
    .value("X_INDX", vpgl_rational_camera<double>::X_INDX)
    .value("Y_INDX", vpgl_rational_camera<double>::Y_INDX)
    .value("Z_INDX", vpgl_rational_camera<double>::Z_INDX)
    .value("U_INDX", vpgl_rational_camera<double>::U_INDX)
    .value("V_INDX", vpgl_rational_camera<double>::V_INDX);

   rational_camera
     .def(py::init<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>,
                   double, double, double, double, double, double, double, double, double, double>())
     .def(py::init<vnl_matrix_fixed<double,4,20>, std::vector<vpgl_scale_offset<double> > >())
     .def("__str__", stream2str<vpgl_rational_camera<double> >)
     .def("coefficient_matrix", &vpgl_rational_camera<double>::coefficient_matrix)
     .def("scale_offsets", &vpgl_rational_camera<double>::scale_offsets)
     .def("project", vpgl_project_point<vpgl_rational_camera<double> >)
     .def("offset", &vpgl_rational_camera<double>::offset)
     .def("project",
      [](vpgl_rational_camera<double> &R, const double x, const double y, const double z)
      {
        vgl_point_3d<double> world_point(x,y,z);                                 
        vgl_point_2d<double> img_point = R.project(world_point);
        return std::make_tuple(img_point.x(), img_point.y());
      }
    );

  m.def("read_rational_camera",
        [](std::string const& fname){return read_rational_camera<double>(fname);},
        py::return_value_policy::take_ownership);

  py::class_<vpgl_local_rational_camera<double>, vpgl_rational_camera<double> >(m, "local_rational_camera")
    .def(py::init<vpgl_lvcs const&, vpgl_rational_camera<double> const&>())
    .def(py::init<double, double, double, vpgl_rational_camera<double> const&>())
    .def("set_lvcs", &vpgl_local_rational_camera<double>::set_lvcs)
    .def("lvcs", &vpgl_local_rational_camera<double>::lvcs)
    .def("project",
      [](vpgl_local_rational_camera<double> &R, const double x, const double y, const double z)
      {
        vgl_point_3d<double> world_point(x,y,z);                                 
        vgl_point_2d<double> img_point = R.project(world_point);
        return std::make_tuple(img_point.x(), img_point.y());
      }
    );

  m.def("read_local_rational_camera",
        [](std::string const& fname){return read_local_rational_camera<double>(fname);},
        py::return_value_policy::take_ownership);


  m.def("compute_affine_from_local_rational", 
     [](vpgl_local_rational_camera<double> &camera, 
      double min_x, double min_y, double min_z, 
      double max_x, double max_y, double max_z, unsigned n_points) 
    {
      if (n_points < 100) {
        n_points = 100; // make it minimum 100 points
      }
      double width = max_x - min_x;
      double depth = max_y - min_y;
      double height = max_z - min_z;
       std::cout << " Using: " << n_points << " points to find the affine camera!\n";
      std::cout << " w: " << width << " d: " << depth << " h: " << height << '\n';
      std::vector< vgl_point_2d<double> > image_pts;
      std::vector< vgl_point_3d<double> > world_pts;
      vnl_random rng;
      for (unsigned i = 0; i < n_points; i++) {
        vgl_point_3d<float> corner_world;
        double x = rng.drand64()*width + min_x;  // sample in local coords
        double y = rng.drand64()*depth + min_y;
        double z = rng.drand64()*height + min_z;
        world_pts.push_back(vgl_point_3d<double>(x,y,z));
        double u, v;
        camera.project(x,y,z,u,v);  // local rational camera has an lvcs, so it handles 
                                    // local coord to global to image point projection internally
        image_pts.push_back(vgl_point_2d<double>(u,v));
      }
      vpgl_affine_camera<double>* out_camera = vpgl_affine_rectification::compute_affine_cam(image_pts, world_pts);
      return out_camera;
    }
  );
  // Geo- Camera definitions
  py::class_<vpgl_geo_camera>(m, "geo_camera")
    // Default methods
    .def(py::init<>())
    .def("__str__", stream2str<vpgl_geo_camera >)
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
    .def("__str__", stream2str<vpgl_lvcs>)

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
            double gx, gy, gz;
            L.local_to_global(lx,ly,lz,ocs,gx,gy,gz,oau,olu);
            return std::make_tuple(gx,gy,gz);
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

}
}
