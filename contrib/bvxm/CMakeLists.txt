cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib-bvxm")

# Add pybind11 module
pybind11_add_module(pybvxm pybvxm.h pybvxm.cxx)

# Link to vxl library
target_link_libraries(pybvxm PRIVATE bvxm)

# Set names
set_target_properties(pybvxm PROPERTIES OUTPUT_NAME "_bvxm")

# install the .so file to the python install dir
install(TARGETS pybvxm DESTINATION ${PYTHON_SITE}/vxl/contrib/bvxm)

# auto generate __init__ file
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/contrib/bvxm/__init__.py \"from ._bvxm import *\")")

