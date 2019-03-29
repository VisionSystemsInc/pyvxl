cmake_minimum_required(VERSION 3.5.1)
project("pyvxl")

# If building PyVXL directly, build VXL
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(PYVXL_TOP_LEVEL TRUE)
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

    # NOTE: If you want to wrap contrib/brl, you'll need to build it in vxl.
    # This currently requires setting the following CMake options:
    # VXL_BUILD_CORE_VIDEO=ON, VXL_BUILD_CONTRIB=ON, VXL_BUILD_BRL=ON

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

# set(PYVXL_LINK_LIBRARIES vgl vnl vpgl vil vgl_algo vpgl_algo vpgl_file_formats)

# Optionally build contrib modules
add_subdirectory("contrib" "pyvxl-contrib-build")

# Add pybind11 modules
pybind11_add_module(pyvnl pyvnl.h pyvnl.cxx)
pybind11_add_module(pyvgl pyvgl.h pyvgl.cxx)
pybind11_add_module(pyvpgl pyvpgl.h pyvpgl.cxx)
pybind11_add_module(pyvil pyvil.h pyvil.cxx)
pybind11_add_module(pyvgl_algo pyvgl_algo.h pyvgl_algo.cxx)
pybind11_add_module(pyvpgl_algo pyvpgl_algo.h pyvpgl_algo.cxx)

# Link
target_link_libraries(pyvnl PRIVATE vnl)
target_link_libraries(pyvgl PRIVATE vgl vnl)
target_link_libraries(pyvpgl PRIVATE vpgl vpgl_file_formats)
target_link_libraries(pyvil PRIVATE vil)
target_link_libraries(pyvgl_algo PRIVATE vgl_algo)
target_link_libraries(pyvpgl_algo PRIVATE vpgl_algo)

# Set names
set_target_properties(pyvnl PROPERTIES OUTPUT_NAME "_vnl")
set_target_properties(pyvgl PROPERTIES OUTPUT_NAME "_vgl")
set_target_properties(pyvpgl PROPERTIES OUTPUT_NAME "_vpgl")
set_target_properties(pyvil PROPERTIES OUTPUT_NAME "_vil")
set_target_properties(pyvgl_algo PROPERTIES OUTPUT_NAME "_vgl_algo")
set_target_properties(pyvpgl_algo PROPERTIES OUTPUT_NAME "_vpgl_algo")

# Find the python install directory
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE_DEFAULT
  OUTPUT_STRIP_TRAILING_WHITESPACE)
set(PYTHON_SITE ${PYTHON_SITE_DEFAULT} CACHE STRING "Python installation directory")

# install the pybind11 extensions
install(TARGETS pyvnl DESTINATION ${PYTHON_SITE}/vxl/vnl)
install(TARGETS pyvgl DESTINATION ${PYTHON_SITE}/vxl/vgl)
install(TARGETS pyvpgl DESTINATION ${PYTHON_SITE}/vxl/vpgl)
install(TARGETS pyvil DESTINATION ${PYTHON_SITE}/vxl/vil)
install(TARGETS pyvgl_algo DESTINATION ${PYTHON_SITE}/vxl/vgl/algo)
install(TARGETS pyvpgl_algo DESTINATION ${PYTHON_SITE}/vxl/vpgl/algo)

# auto generate __init__ files
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vnl/__init__.py \"from ._vnl import *\")")
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vgl/__init__.py \"from ._vgl import *\nfrom .algo._vgl_algo import *\")")
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vpgl/__init__.py \"from ._vpgl import *\nfrom .algo._vpgl_algo import *\")")
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vil/__init__.py \"from ._vil import *\")")
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vgl/algo/__init__.py \"from ._vgl_algo import *\")")
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vpgl/algo/__init__.py \"from ._vpgl_algo import *\")")

# install(TARGETS pyvxl DESTINATION ${PYTHON_SITE})
