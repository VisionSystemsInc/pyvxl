cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-vpgl-file-formats")

# Add pybind11 module
pybind11_add_module(pyvpgl_file_formats pyvpgl_file_formats.h pyvpgl_file_formats.cxx)

# Link to vxl library
target_link_libraries(pyvpgl_file_formats PRIVATE vpgl_file_formats)

# Set names
set_target_properties(pyvpgl_file_formats PROPERTIES OUTPUT_NAME "_vpgl_file_formats")

# install the .so file to the python install dir
install(TARGETS pyvpgl_file_formats DESTINATION ${PYTHON_SITE}/vxl/vpgl/file_formats)

# auto generate __init__ file
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/vpgl/file_formats/__init__.py \"from ._vpgl_file_formats import *\")")

