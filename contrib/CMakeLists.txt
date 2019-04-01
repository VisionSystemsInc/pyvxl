cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib")

# BUILD ALL GUARD
set(PYVXL_CONTRIB_MAKE_ALL FALSE CACHE BOOL "Build every optional vxl contrib module which we have python wrappers for")

# INDIVIDUAL GUARDS, inherit from ALL GUARD
set(PYVXL_CONTRIB_MAKE_BPGL ${PYVXL_CONTRIB_MAKE_ALL} CACHE BOOL "Turn on BPGL build")
set(PYVXL_CONTRIB_MAKE_BRAD ${PYVXL_CONTRIB_MAKE_ALL} CACHE BOOL "Turn on BRAD build")

# Overwrite existing __init__ file
install(CODE "file(WRITE ${PYTHON_SITE}/vxl/contrib/__init__.py \"\")")

if (PYVXL_CONTRIB_MAKE_BPGL)
  # Create the bpgl install directory since it doesn't exist yet
  install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/contrib/bpgl)

  # auto generate __init__ file
  install(CODE "file(APPEND ${PYTHON_SITE}/vxl/contrib/__init__.py \"\nfrom . import bpgl\")")

  # Recurse
  add_subdirectory("bpgl" "bpgl-build")
endif()

if (PYVXL_CONTRIB_MAKE_BRAD)
  # Create the brad install directory since it doesn't exist yet
  install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/contrib/brad)

  # auto generate __init__ file
  install(CODE "file(APPEND ${PYTHON_SITE}/vxl/contrib/__init__.py \"\nfrom . import brad\")")

  # Recurse
  add_subdirectory("brad" "brad-build")
endif()
