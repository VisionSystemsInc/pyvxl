#ifndef pyvil_h_included_
#define pyvil_h_included_

#include <pybind11/pybind11.h>

namespace pyvxl { namespace vil {

void wrap_vil(pybind11::module &m);

}}

#endif
