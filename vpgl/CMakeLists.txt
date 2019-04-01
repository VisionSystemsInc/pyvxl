cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-vpgl")

# Add pybind11 module
pybind11_add_module(pyvpgl pyvpgl.h pyvpgl.cxx)

# Link to vxl library
target_link_libraries(pyvpgl PRIVATE vpgl)

# Set names
set_target_properties(pyvpgl PROPERTIES OUTPUT_NAME "_vpgl")

# install the .so file to the python install dir
install(TARGETS pyvpgl DESTINATION ${PYTHON_SITE}/vxl/vpgl)

# Create the algo and file_formats install directories in case they don't exist
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/vpgl/algo)
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/vpgl/file_formats)

# auto generate __init__ file
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vpgl/__init__.py \"\
from ._vpgl import *\n\
from . import algo\n\
from . import file_formats\
\")")

# Recurse
add_subdirectory("algo" "algo-build")
add_subdirectory("file_formats" "file-formats-build")

