project("pyvxl-vgl")

# add as pyvxl module
pyvxl_add_module(vgl)

# Recurse
add_subdirectory("algo" "algo-build")
