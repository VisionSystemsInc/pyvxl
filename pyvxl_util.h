#ifndef pyvxl_util_h_included
#define pyvxl_util_h_included

#include <string>
#include <sstream>

// Return a string based on output of the object's stream operator
template<typename T>
std::string streamToString(T const& t){

  std::ostringstream buffer;
  buffer << t;
  return buffer.str();
}

#endif
