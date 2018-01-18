project("pyvxl")
include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${PYBIND11_INCLUDE_DIR})

set(PYVXL_SOURCES
  pyvxl.cxx
  pyvnl.h pyvnl.cxx
  pyvgl.h pyvgl.cxx
  pyvpgl.h pyvpgl.cxx
  pyvil.h pyvil.cxx
  pyvgl_algo.h pyvgl_algo.cxx
  )

pybind11_add_module(pyvxl ${PYVXL_SOURCES})
target_link_libraries(pyvxl PRIVATE vnl vil vgl vgl_algo vpgl)
install(TARGETS pyvxl DESTINATION ${PYTHON_LIB_INSTALL_ROOT})
