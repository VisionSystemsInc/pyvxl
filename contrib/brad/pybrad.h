#ifndef pybrad_h_included_
#define pybrad_h_included_

#include <pybind11/pybind11.h>

namespace pyvxl { namespace brad {

void wrap_brad(pybind11::module &m);

}}

#endif
