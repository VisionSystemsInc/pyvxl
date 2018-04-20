project("pyvxl")

include_directories(${PYTHON_INCLUDE_PATH})

set(PYVXL_SOURCES
  pyvxl.cxx
  pyvnl.h pyvnl.cxx
  pyvgl.h pyvgl.cxx
  pyvpgl.h pyvpgl.cxx
  pyvil.h pyvil.cxx
  pyvgl_algo.h pyvgl_algo.cxx
  )

pybind11_add_module(vxl ${PYVXL_SOURCES})
target_link_libraries(vxl PRIVATE vgl vnl vpgl vil)

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE_DEFAULT
  OUTPUT_STRIP_TRAILING_WHITESPACE)
set(PYTHON_SITE ${PYTHON_SITE_DEFAULT} CACHE STRING "Python installation directory")

install(TARGETS vxl DESTINATION ${PYTHON_SITE})
