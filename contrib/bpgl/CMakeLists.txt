cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib-bpgl")

# Create the algo install directory in case it doesn't exist yet
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/contrib/bpgl/algo)

# copy __init__ file over into module
install(FILES __init__.py DESTINATION ${PYTHON_SITE}/vxl/contrib/bpgl)

# Recurse
add_subdirectory("algo" "algo-build")

