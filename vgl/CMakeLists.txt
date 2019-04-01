cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-vgl")

# Add pybind11 module
pybind11_add_module(pyvgl pyvgl.h pyvgl.cxx)

# Link to vxl library
target_link_libraries(pyvgl PRIVATE vgl)

# Set names
set_target_properties(pyvgl PROPERTIES OUTPUT_NAME "_vgl")

# install the .so file to the python install dir
install(TARGETS pyvgl DESTINATION ${PYTHON_SITE}/vxl/vgl)

# Create the algo install directory in case it doesn't exist
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/vgl/algo)

# auto generate __init__ file
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vgl/__init__.py \"from ._vgl import *\nfrom . import algo\")")

# Recurse
add_subdirectory("algo" "algo-build")

