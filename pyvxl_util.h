#ifndef pyvxl_util_h_included
#define pyvxl_util_h_included

#include <string>
#include <sstream>

// Return a string based on output of the object's stream operator
template <class T>
std::string stream2str(T const& obj)
{
  std::stringstream ss;
  ss << obj;
  return ss.str();
}



#endif
