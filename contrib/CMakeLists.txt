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

# is ANY contrib guard active
if (   ${PYVXL_CONTRIB_MAKE_ALL}
    OR ${PYVXL_CONTRIB_MAKE_BPGL}
    OR ${PYVXL_CONTRIB_MAKE_BRAD}
    OR ${PYVXL_CONTRIB_MAKE_BRIP}
    OR ${PYVXL_CONTRIB_MAKE_BVXM}
    OR ${PYVXL_CONTRIB_MAKE_SDET}
    )
  set(pyvxl_contrib_make_any TRUE)
else()
  set(pyvxl_contrib_make_any FALSE)
endif()

# add "contrib" module if necessary
if (pyvxl_contrib_make_any)
  pyvxl_add_module(contrib)

  # Recurse
  if (PYVXL_CONTRIB_MAKE_BPGL)
    add_subdirectory("bpgl" "bpgl-build")
  endif()

  if (PYVXL_CONTRIB_MAKE_BRAD)
    add_subdirectory("brad" "brad-build")
  endif()

  if (PYVXL_CONTRIB_MAKE_BRIP)
    add_subdirectory("brip" "brip-build")
  endif()

  if (PYVXL_CONTRIB_MAKE_BVXM)
    add_subdirectory("bvxm" "bvxm-build")
  endif()

  if (PYVXL_CONTRIB_MAKE_SDET)
    add_subdirectory("sdet" "sdet-build")
  endif()

endif()

