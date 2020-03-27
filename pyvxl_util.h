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


// Conversion from struct to dict
template<typename T>
pybind11::dict struct_to_dict(const T& obj)
{
  pybind11::object pyobj = pybind11::cast(&obj);
  pybind11::dict classDict = pyobj.attr("__class__").attr("__dict__");
  pybind11::dict output;

  for (auto item : classDict)
  {
    switch (PyObject_IsInstance(item.second.ptr(), (PyObject*)(&PyProperty_Type)))
    {
      case 1:
        output[item.first] = pyobj.attr(item.first);
        break;
      case -1:
        throw std::runtime_error("Could not determine the type of "
                                 + item.first.cast<std::string>());
        break;
    }
  }

  return output;
}

// __repr__ via struct_to_dict
template<typename T>
pybind11::str repr_by_dict(const T& obj)
{
  return struct_to_dict(obj).attr("__repr__")();
}

// set object attributes from nested dictionary
void _set_attrs_from_dict(pybind11::object &obj, pybind11::dict dict)
{
  for (auto item : dict) {
    if (pybind11::isinstance<pybind11::dict>(item.second)) {
      auto attr_obj = obj.attr(item.first).cast<pybind11::object>();
      auto value_as_dict = item.second.cast<pybind11::dict>();
      _set_attrs_from_dict(attr_obj, value_as_dict);
    }
    else {
      obj.attr(item.first) = item.second;
    }
  }
}

// set struct values via keyword arguments
template<typename T>
void set_struct_from_kwargs(T &self, pybind11::kwargs kwargs)
{
  pybind11::object pyself = pybind11::cast(&self);
  _set_attrs_from_dict(pyself, kwargs);
}

// init struct via keyword arguments
template<typename T>
T init_struct_from_kwargs(pybind11::kwargs kwargs)
{
  T self;
  set_struct_from_kwargs<T>(self, kwargs);
  return self;
}

#endif
