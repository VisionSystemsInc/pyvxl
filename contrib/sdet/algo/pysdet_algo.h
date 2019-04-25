#ifndef pysdet_algo_h_included_
#define pysdet_algo_h_included_

#include <pybind11/pybind11.h>

namespace pyvxl { namespace sdet { namespace algo {

void wrap_sdet_algo(pybind11::module &m);

}}}

#endif  // pysdet_algo_h_included_
