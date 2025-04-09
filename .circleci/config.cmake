cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)  # same as VXL
project(pyvxl-circleci)

# dependencies
find_package(PythonLibs 3 REQUIRED)

find_package(pybind11 REQUIRED)

# source/binary directories
set(VXL_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/vxl)
set(VXL_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/vxl)

set(PYVXL_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/pyvxl)
set(PYVXL_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/pyvxl)

# VXL build variables
set(VXL_BUILD_CONTRIB ON CACHE BOOL "build vxl contrib")
set(VXL_BUILD_BRL ON CACHE BOOL "build vxl contrib/brl")
set(VXL_BUILD_CORE_VIDEO ON CACHE BOOL "build vxl core/vidl")  # required for brl
set(VXL_BUILD_TESTING OFF CACHE BOOL "build vxl tests")

# pyvxl build variables
set(VXL_DIR ${VXL_SOURCE_DIR})

# add VXL, but only build necessary stuff
add_subdirectory(${VXL_SOURCE_DIR} ${VXL_BINARY_DIR} EXCLUDE_FROM_ALL)

# include core vxl in include paths
set(VXL_CORE_INCLUDE_DIR ${VXL_BINARY_DIR}/core ${VXL_SOURCE_DIR}/core)
set(VXL_VCL_INCLUDE_DIR ${VXL_BINARY_DIR}/vcl ${VXL_SOURCE_DIR}/vcl)

include_directories(${VXL_CORE_INCLUDE_DIR})
include_directories(${VXL_VCL_INCLUDE_DIR})
include_directories(${VXL_DIR}/contrib/brl/bbas)
include_directories(${VXL_DIR}/contrib/brl/bseg)
include_directories(${VXL_DIR}/contrib/gel)

# add pyvxl
add_subdirectory(${PYVXL_SOURCE_DIR} ${PYVXL_BINARY_DIR})
