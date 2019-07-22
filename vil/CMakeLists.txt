cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-vil")

# Add pybind11 module
pybind11_add_module(pyvil pyvil.h pyvil.cxx)

# Link to vxl library
target_link_libraries(pyvil PRIVATE vil)

# Set names
set_target_properties(pyvil PROPERTIES OUTPUT_NAME "_vil")

# install the .so file to the python install dir
install(TARGETS pyvil DESTINATION ${PYTHON_SITE}/vxl/vil)

# copy __init__ file over into module
install(FILES __init__.py DESTINATION ${PYTHON_SITE}/vxl/vil)

