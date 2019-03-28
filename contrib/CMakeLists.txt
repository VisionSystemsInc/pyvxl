cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib")

# BUILD ALL GUARD
set(PYVXL_CONTRIB_MAKE_ALL FALSE CACHE BOOL "Build every optional vxl contrib module which we have python wrappers for")

# INDIVIDUAL GUARDS, inherit from ALL GUARD
set(PYVXL_CONTRIB_MAKE_BPGL ${PYVXL_CONTRIB_MAKE_ALL} CACHE BOOL "Turn on BPGL build")

# Initialize sources and link libs as empty
set(PYVXL_CONTRIB_SOURCES pycontrib.cxx)
set(PYVXL_CONTRIB_LINK_LIBRARIES)


message(STATUS "in contrib cmakelist!")
if (PYVXL_CONTRIB_MAKE_BPGL)
  message(STATUS "pyvxl contrib: wrapping contrib/bpgl libraries")
  add_definitions("-DPYVXL_WITH_CONTRIB_BPGL")
  set(PYVXL_CONTRIB_SOURCES ${PYVXL_CONTRIB_SOURCES} pybpgl_algo.h pybpgl_algo.cxx)
  set(PYVXL_CONTRIB_LINK_LIBRARIES ${PYVXL_CONTRIB_LINK_LIBRARIES} bpgl_algo)
  if (PYVXL_TOP_LEVEL)
    include_directories(${VXL_DIR}/contrib/brl/bbas)
  endif()
endif()

pybind11_add_module(pyvxl-contrib ${PYVXL_CONTRIB_SOURCES})

target_link_libraries(pyvxl-contrib PRIVATE ${PYVXL_CONTRIB_LINK_LIBRARIES})
set_target_properties(pyvxl-contrib PROPERTIES OUTPUT_NAME "_vxl_contrib")

install(TARGETS pyvxl-contrib DESTINATION ${PYTHON_SITE})

