cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-vgl-algo")

# Add pybind11 module
pybind11_add_module(pyvgl_algo pyvgl_algo.h pyvgl_algo.cxx)

# Link to vxl library
target_link_libraries(pyvgl_algo PRIVATE vgl_algo)

# Set names
set_target_properties(pyvgl_algo PROPERTIES OUTPUT_NAME "_vgl_algo")

# install the .so file to the python install dir
install(TARGETS pyvgl_algo DESTINATION ${PYTHON_SITE}/vxl/vgl/algo)

# auto generate __init__ file
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vgl/algo/__init__.py \"from ._vgl_algo import *\")")

