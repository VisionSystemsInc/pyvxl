cmake_minimum_required(VERSION 3.5.1)
project("pyvxl")

# pyvxl project source directory
set(PYVXL_PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}")

# print pyvxl debug messages during cmake
set(PYVXL_DEBUG_MESSAGES "OFF" CACHE BOOL "Print pyvxl debug messages")

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
           "${PYTHON_SITE}/vxl"          # replace with this
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
    if (EXISTS "${target_file_h}")
        message(STATUS "pyvxl found ${target_name}")

        # Add pybind11 module
        pybind11_add_module(${target_name} ${target_files})

        # Link to vxl library
        target_link_libraries(${target_name} PRIVATE ${vxl_name})

        # Set output_name of target
        set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${output_name})

        # target installation
        install(TARGETS ${target_name} DESTINATION ${install_dir})

    endif()

endfunction()


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

include_directories(${PROJECT_SOURCE_DIR})

# Find the python install directory
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE_DEFAULT
  OUTPUT_STRIP_TRAILING_WHITESPACE)
set(PYTHON_SITE ${PYTHON_SITE_DEFAULT} CACHE STRING "Python installation directory")

# add main module
pyvxl_add_module(vxl)

# Recurse
add_subdirectory("vgl" "pyvxl_build/vgl-build")
add_subdirectory("vil" "pyvxl_build/vil-build")
add_subdirectory("vnl" "pyvxl_build/vnl-build")
add_subdirectory("vpgl" "pyvxl_build/vpgl-build")
add_subdirectory("contrib" "pyvxl_build/contrib-build")

