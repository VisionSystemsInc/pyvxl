#ifndef acal_pystruct_h_included_
#define acal_pystruct_h_included_

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace py::literals;


namespace pyvxl { namespace acal {

template<typename T>
  T dict_to_struct(const py::dict& d){

    T obj;
    py::object pyobj = py::cast(&obj);

    for(auto item : d){
      pyobj.attr(item.first) = item.second;
    }

    return obj;
  }

  template<typename T>
  py::dict struct_to_dict(const T& obj){

    py::object pyobj = py::cast(&obj);
    py::dict classDict = pyobj.attr("__class__").attr("__dict__");
    py::dict output;

    for(auto item : classDict){
      switch(PyObject_IsInstance(item.second.ptr(), (PyObject*)(&PyProperty_Type))){
      case 1:
        output[item.first] = pyobj.attr(item.first);
        break;
      case -1:
        throw std::runtime_error("Could not determine the type of " + item.first.cast<std::string>());
        break;
      }
    }

    return output;
  }

  template<typename T>
  py::str repr_by_dict(const T& obj) {
    return struct_to_dict(obj).attr("__repr__")();
  }

}}

#endif
