#ifndef pybsgm_h_included_
#define pybsgm_h_included_

#include <pybind11/pybind11.h>

namespace pyvxl { namespace bsgm {

void wrap_bsgm(pybind11::module &m);

}}

#endif  // pybsgm_h_included_
