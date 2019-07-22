cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-vpgl")

# Add pybind11 module
pybind11_add_module(pyvpgl pyvpgl.h pyvpgl.cxx)

# Link to vxl library
target_link_libraries(pyvpgl PRIVATE vpgl vpgl_file_formats brip vsol)

# Set names
set_target_properties(pyvpgl PROPERTIES OUTPUT_NAME "_vpgl")

# install the .so file to the python install dir
install(TARGETS pyvpgl DESTINATION ${PYTHON_SITE}/vxl/vpgl)

# Create the algo install directory in case it doesn't exist
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/vpgl/algo)

# copy __init__ file over into module
install(FILES __init__.py DESTINATION ${PYTHON_SITE}/vxl/vpgl)

# Recurse
add_subdirectory("algo" "algo-build")

