project("pyvxl-vbl")

# add as pyvxl module
pyvxl_add_module(vbl)

# additional linking
target_link_libraries(pyvbl PRIVATE vbl_io)
