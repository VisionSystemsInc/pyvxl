cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-vgl")

# add as pyvxl module
pyvxl_add_module(vgl)

# Recurse
add_subdirectory("algo" "algo-build")
