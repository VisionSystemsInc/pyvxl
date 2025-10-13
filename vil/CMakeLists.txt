project("pyvxl-vil")

# add as pyvxl module
pyvxl_add_module(vil)

# additional linking
target_link_libraries(pyvil PRIVATE vil_io)
