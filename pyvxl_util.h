#ifndef pyvxl_util_h_included
#define pyvxl_util_h_included

#include <pybind11/pybind11.h>
#include <string>
#include <sstream>
#include <vsl/vsl_binary_io.h>

// Return a string based on output of the object's stream operator
template<typename T>
std::string streamToString(T const& t){

  std::ostringstream buffer;
  buffer << t;
  return buffer.str();
}

/* PICKLE VXL CORE INSTANCE with vsl_binary_io
 *
 * Any VXL core class with an "io" counterpart can be pickled via
 * these helper functions.
 *
 * Pickle capabilites are added to an object via py::pickle.
 * Note the appropriate "io" header must also be made available.
 *
 *   -----
 *   #include "../pyvxl_util.h"
 *   #include <vpgl/io/vpgl_io_affine_camera.h>
 *
 *   py::class_< vpgl_affine_camera<double> > (m, "affine_camera")
 *     .def(py::init ...)
 *     ...
 *     .def(py::pickle(&vslPickleGetState< vpgl_affine_camera<double> >,
 *                     &vslPickleSetState< vpgl_affine_camera<double> >))
 *     ;
 *   -----
 *
 * Additionally, "target_link_libraries" in CMakeLists.txt must link
 * to the appropriate "io" target
 *
 *   -----
 *   target_link_libraries(pyvpgl PRIVATE vpgl vpgl_io ...)
 *                                             ^^^^^^^
 *   -----
 *
 * Assuming the availability of "operator==" for the pickled class,
 * a python pickle unit test follows this pattern:
 *
 *   -----
 *   def test_pickle(self):
 *     objA = self.init_obj()
 *     objB = pickle.loads(pickle.dumps(objA))
 *     self.assertEqual(objA, objB)
 *   -----
 */

template<typename T>
pybind11::bytes vslPickleGetState(T const& obj)
{
  std::ostringstream oss;
  vsl_b_ostream oss_vsl(&oss);
  vsl_b_write(oss_vsl, obj);
  return pybind11::bytes(oss.str());
}

template<typename T>
T vslPickleSetState(pybind11::bytes b)
{
  T obj;
  std::istringstream iss(b);
  vsl_b_istream iss_vsl(&iss);
  vsl_b_read(iss_vsl, obj);
  return obj;
}

#endif
