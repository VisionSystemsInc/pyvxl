cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-vpgl-algo")

# Add pybind11 module
pybind11_add_module(pyvpgl_algo pyvpgl_algo.h pyvpgl_algo.cxx)

# Link to vxl library
target_link_libraries(pyvpgl_algo PRIVATE vpgl_algo)

# Set names
set_target_properties(pyvpgl_algo PROPERTIES OUTPUT_NAME "_vpgl_algo")

# install the .so file to the python install dir
install(TARGETS pyvpgl_algo DESTINATION ${PYTHON_SITE}/vxl/vpgl/algo)

# auto generate __init__ file
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vpgl/algo/__init__.py \"from ._vpgl_algo import *\")")

