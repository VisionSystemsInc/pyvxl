cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-vpgl")

# add as pyvxl module
pyvxl_add_module(vpgl)

# additional linking
target_link_libraries(pyvpgl PRIVATE vpgl_io vpgl_file_formats)

# Recurse
add_subdirectory("algo" "algo-build")
