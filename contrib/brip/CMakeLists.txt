cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib-brip")

# Add pybind11 module
pybind11_add_module(pybrip pybrip.h pybrip.cxx)

# Link to vxl library
target_link_libraries(pybrip PRIVATE brip)

# Set names
set_target_properties(pybrip PROPERTIES OUTPUT_NAME "_brip")

# install the .so file to the python install dir
install(TARGETS pybrip DESTINATION ${PYTHON_SITE}/vxl/contrib/brip)

# copy __init__ file over into module
install(FILES __init__.py DESTINATION ${PYTHON_SITE}/vxl/contrib/brip)

