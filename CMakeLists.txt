include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${PYBIND11_INCLUDE_PATH})

set(PYVXL_SOURCES
  pyvxl.cxx
  pyvnl.h pyvnl.cxx
  pyvgl.h pyvgl.cxx
  pyvpgl.h pyvpgl.cxx
  pyvil.h pyvil.cxx
  pyvgl_algo.h pyvgl_algo.cxx
  )

# pybind11_add_module
add_library(pyvxl SHARED ${PYVXL_SOURCES})
target_link_libraries(pyvxl vgl vnl vpgl vil ${PYTHON_LIBRARIES})
set_target_properties(pyvxl PROPERTIES OUTPUT_NAME vxl)
set_target_properties(pyvxl PROPERTIES PREFIX "")
set_target_properties(pyvxl PROPERTIES SUFFIX ".so")


install(TARGETS pyvxl DESTINATION ${PYTHON_LIB_INSTALL_ROOT})
