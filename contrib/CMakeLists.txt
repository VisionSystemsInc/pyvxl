# contrib modules are entirely optional, controlled by boolean cache
# guards of the form PYVXL_CONTRIB_MAKE_*
#
# PYVXL_CONTRIB_MAKE_ALL sets the default across all contrib directories
#
# Users may also specify per-module inclusion, for example:
#   PYVXL_CONTRIB_MAKE_BPGL controls BPGL inclusion
#   PYVXL_CONTRIB_MAKE_SDET controls SDET inclusion
#   etc.
#
# For developers, this CMakeLists.txt file will auto-discover any available
# contrib modules, initialze a PYVXL_CONTRIB_MAKE_* guard if one does not
# exist, and add the contrib subdirectory if necessary.
#
project("pyvxl-contrib")

# build all guard
set(PYVXL_CONTRIB_MAKE_ALL FALSE CACHE BOOL "Build every optional pyvxl contrib module")

# list of contrib components
set(contrib_components bpgl brad brip bvxm sdet)

# has base contrib module already been added?
set(contrib_added FALSE)

# add contrib components
foreach(component ${contrib_components})

  # cache variable for boolean guard
  string(TOUPPER "${component}" component_upper)
  set(guard_name "PYVXL_CONTRIB_MAKE_${component_upper}")
  set(${guard_name} ${PYVXL_CONTRIB_MAKE_ALL} CACHE BOOL "Turn on ${component_upper} build")

  # add component if requested
  if(${${guard_name}})

    # add base contrib directory
    if(NOT contrib_added)
      pyvxl_add_module(contrib)
      set(contrib_added TRUE)
    endif()

    # add component subdirectory
    add_subdirectory("${component}" "${component}-build")

  endif()

endforeach()
