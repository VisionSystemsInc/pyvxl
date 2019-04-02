#ifndef pyvnl_h_included_
#define pyvnl_h_included_

#include <pybind11/pybind11.h>

namespace pyvxl { namespace vnl {

void wrap_vnl(pybind11::module &m);

}}

#endif
