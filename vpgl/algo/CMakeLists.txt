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

# copy __init__ file over into module
install(FILES __init__.py DESTINATION ${PYTHON_SITE}/vxl/vpgl/algo)

