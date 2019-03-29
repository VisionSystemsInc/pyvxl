cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib")

# BUILD ALL GUARD
set(PYVXL_CONTRIB_MAKE_ALL FALSE CACHE BOOL "Build every optional vxl contrib module which we have python wrappers for")

# INDIVIDUAL GUARDS, inherit from ALL GUARD
set(PYVXL_CONTRIB_MAKE_BPGL ${PYVXL_CONTRIB_MAKE_ALL} CACHE BOOL "Turn on BPGL build")

# Initialize sources and link libs as empty
# set(PYVXL_CONTRIB_SOURCES pycontrib.cxx)
# set(PYVXL_CONTRIB_LINK_LIBRARIES)

# Find the python install directory
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE_DEFAULT
  OUTPUT_STRIP_TRAILING_WHITESPACE)
set(PYTHON_SITE ${PYTHON_SITE_DEFAULT} CACHE STRING "Python installation directory")

install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/contrib)
# auto generate __init__ file
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/contrib/__init__.py \"\")")

if (PYVXL_CONTRIB_MAKE_BPGL)
  message(STATUS "pyvxl contrib: wrapping contrib/bpgl libraries")
  add_definitions("-DPYVXL_WITH_CONTRIB_BPGL")

  # Add pybind11 module
  pybind11_add_module(pybpgl_algo pybpgl_algo.h pybpgl_algo.cxx)

  # Link
  target_link_libraries(pybpgl_algo PRIVATE bpgl_algo)
  # set(PYVXL_CONTRIB_LINK_LIBRARIES ${PYVXL_CONTRIB_LINK_LIBRARIES} bpgl_algo)

  # Set name
  set_target_properties(pybpgl_algo PROPERTIES OUTPUT_NAME "_bpgl_algo")

  # install the pybind11 extension
  install(TARGETS pybpgl_algo DESTINATION ${PYTHON_SITE}/vxl/contrib/bpgl/algo)

  # append to base __init__
  install(CODE "file(APPEND ${PYTHON_SITE}/vxl/contrib/__init__.py \"\nfrom . import bpgl\")")

  # auto generate __init__ file
  install(CODE "file(WRITE ${PYTHON_SITE}/vxl/contrib/bpgl/__init__.py \"from . import algo\")")
  install(CODE "file(WRITE ${PYTHON_SITE}/vxl/contrib/bpgl/algo/__init__.py \"from ._bpgl_algo import *\")")

  if (PYVXL_TOP_LEVEL)
    include_directories(${VXL_DIR}/contrib/brl/bbas)
  endif()
endif()

