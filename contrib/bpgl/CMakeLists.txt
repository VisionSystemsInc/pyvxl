cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib-bpgl")

# Create the algo install directory in case it doesn't exist yet
install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/contrib/bpgl/algo)

# auto generate __init__ file to load algo submodule by default
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/contrib/bpgl/__init__.py \"from . import algo\")")

# Recurse
add_subdirectory("algo" "algo-build")

