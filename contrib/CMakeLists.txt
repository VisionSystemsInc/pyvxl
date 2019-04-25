cmake_minimum_required(VERSION 3.5.1)
project("pyvxl-contrib")

# BUILD ALL GUARD
set(PYVXL_CONTRIB_MAKE_ALL FALSE CACHE BOOL "Build every optional vxl contrib module which we have python wrappers for")

# INDIVIDUAL GUARDS, inherit from ALL GUARD
set(PYVXL_CONTRIB_MAKE_BPGL ${PYVXL_CONTRIB_MAKE_ALL} CACHE BOOL "Turn on BPGL build")
set(PYVXL_CONTRIB_MAKE_BRAD ${PYVXL_CONTRIB_MAKE_ALL} CACHE BOOL "Turn on BRAD build")
set(PYVXL_CONTRIB_MAKE_BRIP ${PYVXL_CONTRIB_MAKE_ALL} CACHE BOOL "Turn on BRIP build")
set(PYVXL_CONTRIB_MAKE_BVXM ${PYVXL_CONTRIB_MAKE_ALL} CACHE BOOL "Turn on BVXM build")
set(PYVXL_CONTRIB_MAKE_SDET ${PYVXL_CONTRIB_MAKE_ALL} CACHE BOOL "Turn on SDET build")

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

if (PYVXL_CONTRIB_MAKE_BRIP)
  # Create the brip install directory since it doesn't exist yet
  install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/contrib/brip)

  # auto generate __init__ file
  install(CODE "file(APPEND ${PYTHON_SITE}/vxl/contrib/__init__.py \"\nfrom . import brip\")")

  # Recurse
  add_subdirectory("brip" "brip-build")
endif()

if (PYVXL_CONTRIB_MAKE_BVXM)
  # Create the bvxm install directory since it doesn't exist yet
  install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/contrib/bvxm)

  # auto generate __init__ file
  install(CODE "file(APPEND ${PYTHON_SITE}/vxl/contrib/__init__.py \"\nfrom . import bvxm\")")

  # Recurse
  add_subdirectory("bvxm" "bvxm-build")
endif()

if (PYVXL_CONTRIB_MAKE_SDET)
  # Create the sdet install directory since it doesn't exist yet
  install(DIRECTORY DESTINATION ${PYTHON_SITE}/vxl/contrib/sdet)

  # auto generate __init__ file
  install(CODE "file(APPEND ${PYTHON_SITE}/vxl/contrib/__init__.py \"\nfrom . import sdet\")")

  # Recurse
  add_subdirectory("sdet" "sdet-build")
endif()

