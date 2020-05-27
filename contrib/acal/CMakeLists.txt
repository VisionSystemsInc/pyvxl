project("pyvxl-contrib-acal")

# add as pyvxl module
pyvxl_add_module(acal)

# additional linking
target_link_libraries(pyacal PRIVATE acal_io)
