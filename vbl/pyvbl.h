#ifndef pyvbl_h_included_
#define pyvbl_h_included_

#include <pybind11/pybind11.h>

namespace pyvxl { namespace vbl {

void wrap_vbl(pybind11::module &m);

}}

#endif
