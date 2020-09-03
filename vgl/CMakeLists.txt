project("pyvxl-vgl")

# add as pyvxl module
pyvxl_add_module(vgl)

# additional linking
target_link_libraries(pyvgl PRIVATE vgl_io)

# Recurse
add_subdirectory("algo" "algo-build")
