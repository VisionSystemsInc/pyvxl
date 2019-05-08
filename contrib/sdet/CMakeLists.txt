cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib-sdet")

# Add pybind11 module
pybind11_add_module(pysdet pysdet.h pysdet.cxx)

# Link to vxl library
target_link_libraries(pysdet PRIVATE sdet)

# Set names
set_target_properties(pysdet PROPERTIES OUTPUT_NAME "_sdet")

# install the .so file to the python install dir
install(TARGETS pysdet DESTINATION ${PYTHON_SITE}/vxl/contrib/sdet)

# copy __init__ file over into module
install(FILES __init__.py DESTINATION ${PYTHON_SITE}/vxl/contrib/sdet)

# Recurse
add_subdirectory("algo" "algo-build")

