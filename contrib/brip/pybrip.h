#ifndef pybrip_h_included_
#define pybrip_h_included_

#include <pybind11/pybind11.h>

namespace pyvxl { namespace brip {

void wrap_brip(pybind11::module &m);

}}

#endif  // pybrip_h_included_
