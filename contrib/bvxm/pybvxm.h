#ifndef pybvxm_h_included_
#define pybvxm_h_included_

#include <pybind11/pybind11.h>

namespace pyvxl { namespace bvxm {

void wrap_bvxm(pybind11::module &m);

}}

#endif  // pybvxm_h_included_
