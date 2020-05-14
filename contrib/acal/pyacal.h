#ifndef pyacal_h_included_
#define pyacal_h_included_

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace pyvxl { namespace acal {

void wrap_f_params(py::module &m);
void wrap_match_params(py::module &m);
void wrap_acal_corr(py::module &m);
void wrap_acal_match_pair(py::module &m);
void wrap_match_vertex(py::module &m);
void wrap_match_edge(py::module &m);
void wrap_acal_match_tree(py::module &m);
void wrap_acal_match_graph(py::module &m);

}}

#endif  // pyacal_h_included_
