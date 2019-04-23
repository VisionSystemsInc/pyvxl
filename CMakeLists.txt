cmake_minimum_required(VERSION 3.5.1)
project("pyvxl")

# If building PyVXL directly, build VXL and pybind11
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
    # set(VXL_BUILD_CONTRIB "TRUE" CACHE BOOL "BUILD CONTRIB")
    # set(VXL_BUILD_CORE_VIDEO "TRUE" CACHE BOOL "BUILD CORE VIDEO")
    # set(VXL_BUILD_BRL "TRUE" CACHE BOOL "BUILD BRL")

    # Build VXL
    add_subdirectory(${VXL_DIR} ${VXL_BINARY_DIR} EXCLUDE_FROM_ALL)

    # include vxl headers
    include_directories(${VXL_DIR})
    include_directories(${VXL_DIR}/core ${VXL_BINARY_DIR}/core)
    include_directories(${VXL_DIR}/vcl ${VXL_BINARY_DIR}/vcl)
    include_directories(${VXL_BINARY_DIR}/contrib)
    include_directories(${VXL_DIR}/contrib/brl)
    include_directories(${VXL_DIR}/contrib/brl/bseg)
    include_directories(${VXL_DIR}/contrib/brl/bbas)
    include_directories(${VXL_DIR}/contrib/gel)

    # Build Pybind11
    set(PYBIND11_DIR "PYBIND11_NOTFOUND" CACHE PATH "Location of Pybind11 source directory")
    message(STATUS "PYBIND11_DIR = ${PYBIND11_DIR}")
    if (NOT IS_DIRECTORY ${PYBIND11_DIR})
      message(FATAL_ERROR "Set PYBIND11_DIR to the location of the Pybind11 source directory")
    endif()
    add_subdirectory(${PYBIND11_DIR} "pybind_build")
endif()

include_directories(${PYTHON_INCLUDE_PATH})

# Find the python install directory
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE_DEFAULT
  OUTPUT_STRIP_TRAILING_WHITESPACE)
set(PYTHON_SITE ${PYTHON_SITE_DEFAULT} CACHE STRING "Python installation directory")

# Create the install directories in case they aren't installed to, so auto-loading doesn't break
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/vgl)
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/vil)
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/vnl)
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/vpgl)
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/contrib)

# generate __init__ file to auto-load modules by default
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/__init__.py \"\
from . import vgl\n\
from . import vil\n\
from . import vnl\n\
from . import vpgl\n\
from . import contrib\n\
\")")

# Recurse
add_subdirectory("vgl" "vgl-build")
add_subdirectory("vil" "vil-build")
add_subdirectory("vnl" "vnl-build")
add_subdirectory("vpgl" "vpgl-build")
add_subdirectory("contrib" "contrib-build")

