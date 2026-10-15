cmake_minimum_required(VERSION 3.21)
project("pyvxl")

if(NOT CMAKE_CXX_STANDARD)
   set(CMAKE_CXX_STANDARD "17")
endif()

# pyvxl project source directory
set(PYVXL_PROJECT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# print pyvxl debug messages during cmake
set(PYVXL_DEBUG_MESSAGES "OFF" CACHE BOOL "Print pyvxl debug messages")

# header-only target for top-level utilities
add_library(pyvxl_utils INTERFACE)
target_include_directories(pyvxl_utils INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
add_library(pyvxl::utils ALIAS pyvxl_utils)

# function: add pyvxl module
function(pyvxl_add_module vxl_name)

    # target variables
    set(target_name "py${vxl_name}")
    set(output_name "_${vxl_name}")

    # target files
    set(target_file_h "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.h")
    set(target_file_cxx "${CMAKE_CURRENT_SOURCE_DIR}/${target_name}.cxx")
    set(target_files ${target_file_h} ${target_file_cxx})

    # install directory
    string(REPLACE
           "${PYVXL_PROJECT_SOURCE_DIR}" # string to be replaced
           "${PYVXL_INSTALL_DIR}"          # replace with this
           install_dir                   # output
           "${CMAKE_CURRENT_SOURCE_DIR}" # input
           )

    # debug messages (currently disabled)
    if (PYVXL_DEBUG_MESSAGES)
        message(STATUS "PYVXL_ADD_MODULE:")
        message(STATUS "  vxl_name     = ${vxl_name}")
        message(STATUS "  target_name  = ${target_name}")
        message(STATUS "  output_name  = ${output_name}")
        message(STATUS "  target_files = ${target_files}")
        message(STATUS "  install_dir  = ${install_dir}")
    endif()

    # install directory setup
    install(DIRECTORY DESTINATION ${install_dir})
    install(FILES __init__.py DESTINATION ${install_dir})

    # initialize target if necessary
    if (EXISTS "${target_file_cxx}")
        message(STATUS "pyvxl found ${target_name}")

        # Add pybind11 module
        pybind11_add_module(${target_name} ${target_files})

        # Link to vxl library & utilities
        target_link_libraries(${target_name} PRIVATE ${vxl_name} pyvxl::utils)

        # Set output_name of target
        set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${output_name})

        # target installation
        install(TARGETS ${target_name} DESTINATION ${install_dir})

    endif()

endfunction()

# python
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# pybind11 via find_package or add_subdirectory
# https://pybind11.readthedocs.io/en/stable/compiling.html#find-package-vs-add-subdirectory
set(PYBIND11_SOURCE_DIR "PYBIND11_SOURCE_DIR-NOTFOUND" CACHE PATH "Location of pybind11 source directory")
if (PYBIND11_SOURCE_DIR)
  message(STATUS "pybind11 source directory at ${PYBIND11_SOURCE_DIR}")
  if (NOT IS_DIRECTORY "${PYBIND11_SOURCE_DIR}")
    message(FATAL_ERROR "Set PYBIND11_SOURCE_DIR to the location of the pybind11 source directory")
  endif()
  add_subdirectory("${PYBIND11_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/pybind11_build")
else()
  find_package(pybind11 REQUIRED)
endif()

# If building PyVXL directly, also include VXL
if(PROJECT_IS_TOP_LEVEL)

    # find and include VXL
    set(VXL_SOURCE_DIR "VXL_SOURCE_DIR-NOTFOUND" CACHE PATH "Location of VXL source directory")
    message(STATUS "VXL_SOURCE_DIR = ${VXL_SOURCE_DIR}")
    if (NOT IS_DIRECTORY ${VXL_SOURCE_DIR})
      message(FATAL_ERROR "Set VXL_SOURCE_DIR to the location of VXL source directory")
    endif()
    set(VXL_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/vxl_build")

    # NOTE: If you want to wrap contrib/brl, you'll need to build it in vxl.
    # This currently requires setting the following CMake options:
    # set(VXL_BUILD_CONTRIB "TRUE" CACHE BOOL "BUILD CONTRIB")
    # set(VXL_BUILD_CORE_VIDEO "TRUE" CACHE BOOL "BUILD CORE VIDEO")
    # set(VXL_BUILD_BRL "TRUE" CACHE BOOL "BUILD BRL")

    # Build VXL
    add_subdirectory(${VXL_SOURCE_DIR} ${VXL_BINARY_DIR} EXCLUDE_FROM_ALL)

    # include vxl headers
    include_directories(${VXL_SOURCE_DIR})
    include_directories(${VXL_SOURCE_DIR}/core ${VXL_BINARY_DIR}/core)
    include_directories(${VXL_SOURCE_DIR}/vcl ${VXL_BINARY_DIR}/vcl)
    include_directories(${VXL_BINARY_DIR}/contrib)
    include_directories(${VXL_SOURCE_DIR}/contrib/brl)
    include_directories(${VXL_SOURCE_DIR}/contrib/brl/bseg)
    include_directories(${VXL_SOURCE_DIR}/contrib/brl/bbas)
    include_directories(${VXL_SOURCE_DIR}/contrib/gel)

endif()

# python site-packages cache variable
set(PYTHON_SITE_DEFAULT "${Python_SITELIB}")
set(PYTHON_SITE "${PYTHON_SITE_DEFAULT}" CACHE STRING "Python site-packages directory")

if ("${PYTHON_SITE}" STREQUAL "${PYTHON_SITE_DEFAULT}")
  set(PYTHON_SITE_INITIALIZED_TO_DEFAULT 1)
endif()

if (PYVXL_DEBUG_MESSAGES)
    message(STATUS "PYTHON_SITE_DEFAULT=${PYTHON_SITE_DEFAULT}")
    message(STATUS "PYTHON_SITE=${PYTHON_SITE}")
    message(STATUS "PYTHON_SITE_INITIALIZED_TO_DEFAULT=${PYTHON_SITE_INITIALIZED_TO_DEFAULT}")
endif()

# installation directory as sub-dir of PYTHON_SITE
set(PYVXL_INSTALL_DIR "${PYTHON_SITE}/vxl")
message(STATUS "pyvxl will be installed to ${PYVXL_INSTALL_DIR}")

# add main module
pyvxl_add_module(vxl)

# add core components
set(core_components vbl vgl vil vnl vpgl)
foreach(component ${core_components})
    add_subdirectory("${component}" "pyvxl_build/${component}-build")
endforeach()

# add contrib components
add_subdirectory("contrib" "pyvxl_build/contrib-build")
