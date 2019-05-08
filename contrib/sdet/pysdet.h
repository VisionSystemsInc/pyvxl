#ifndef pysdet_h_included_
#define pysdet_h_included_

#include <pybind11/pybind11.h>

namespace pyvxl { namespace sdet {

void wrap_sdet(pybind11::module &m);

}}

#endif  // pysdet_h_included_
