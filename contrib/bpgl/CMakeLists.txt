cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib-bpgl")

# add as pyvxl module
pyvxl_add_module(bpgl)

# Recurse
add_subdirectory("algo" "algo-build")
