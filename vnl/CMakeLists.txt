cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-vnl")

# Add pybind11 module
pybind11_add_module(pyvnl pyvnl.h pyvnl.cxx)

# Link to vxl library
target_link_libraries(pyvnl PRIVATE vnl)

# Set names
set_target_properties(pyvnl PROPERTIES OUTPUT_NAME "_vnl")

# install the .so file to the python install dir
install(TARGETS pyvnl DESTINATION ${PYTHON_SITE}/vxl/vnl)

# auto generate __init__ file
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vnl/__init__.py \"from ._vnl import *\")")

