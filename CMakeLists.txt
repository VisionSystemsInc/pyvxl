cmake_minimum_required(VERSION 3.5.1)
project("pyvxl")
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    # pyvxl built as top-level project

    # find Python
    find_package(PythonLibs 3 REQUIRED)

    # find and include VXL
    set(VXL_DIR "VXL_NOTFOUND" CACHE PATH "Location of VXL source directory")
    message(STATUS "VXL_DIR = ${VXL_DIR}")
    if (NOT IS_DIRECTORY ${VXL_DIR})
      message(FATAL_ERROR "Set VXL_DIR to the location of VXL source directory")
    endif()
    set(VXL_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/vxl_build)
    add_subdirectory(${VXL_DIR} ${VXL_BINARY_DIR} EXCLUDE_FROM_ALL)

    # include core vxl in include paths
    set(VXL_CORE_INCLUDE_DIR ${VXL_BINARY_DIR}/core ${VXL_DIR}/core)
    set(VXL_VCL_INCLUDE_DIR ${VXL_BINARY_DIR}/vcl ${VXL_DIR}/vcl)
    include_directories(${VXL_CORE_INCLUDE_DIR})
    include_directories(${VXL_VCL_INCLUDE_DIR})

    set(PYBIND11_DIR "PYBIND11_NOTFOUND" CACHE PATH "Location of Pybind11 source directory")
    message(STATUS "PYBIND11_DIR = ${PYBIND11_DIR}")
    if (NOT IS_DIRECTORY ${PYBIND11_DIR})
      message(FATAL_ERROR "Set PYBIND11_DIR to the location of the Pybind11 source directory")
    endif()
    add_subdirectory(${PYBIND11_DIR} "pybind_build")
endif()

include_directories(${PYTHON_INCLUDE_PATH})

set(PYVXL_SOURCES
  pyvxl.cxx
  pyvnl.h pyvnl.cxx
  pyvgl.h pyvgl.cxx
  pyvpgl.h pyvpgl.cxx
  pyvil.h pyvil.cxx
  pyvgl_algo.h pyvgl_algo.cxx
  pyvpgl_algo.h pyvpgl_algo.cxx
  )

pybind11_add_module(vxl ${PYVXL_SOURCES})
target_link_libraries(vxl PRIVATE vgl vnl vpgl vil vgl_algo vpgl_algo)

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE_DEFAULT
  OUTPUT_STRIP_TRAILING_WHITESPACE)
set(PYTHON_SITE ${PYTHON_SITE_DEFAULT} CACHE STRING "Python installation directory")

install(TARGETS vxl DESTINATION ${PYTHON_SITE})
