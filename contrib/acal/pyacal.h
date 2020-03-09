#ifndef pyacal_h_included_
#define pyacal_h_included_

#include <pybind11/pybind11.h>

namespace pyvxl { namespace acal {

void wrap_acal(pybind11::module &m);

}}

#endif  // pyacal_h_included_
