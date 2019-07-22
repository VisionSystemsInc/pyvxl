cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib-bvxm-algo")

# Add pybind11 module
pybind11_add_module(pybvxm_algo pybvxm_algo.h pybvxm_algo.cxx)

# Link to vxl library
target_link_libraries(pybvxm_algo PRIVATE bvxm_algo)

# Set names
set_target_properties(pybvxm_algo PROPERTIES OUTPUT_NAME "_bvxm_algo")

# install the .so file to the python install dir
install(TARGETS pybvxm_algo DESTINATION ${PYTHON_SITE}/vxl/contrib/bvxm/algo)

# copy __init__ file over into module
install(FILES __init__.py DESTINATION ${PYTHON_SITE}/vxl/contrib/bvxm/algo)

