# pyvxl
Python wrapper for vxl using pybind11

[![CircleCI](https://circleci.com/gh/VisionSystemsInc/pyvxl.svg?style=svg)](https://circleci.com/gh/VisionSystemsInc/pyvxl)

# Requirements

- pybind11 - You need to add the pybind11 repo as a submodule to your project (or install it into the OS)
- vxl

# Cmake

pyvxl can be build standalone or as a module of an existing project.

## Contrib
Add 

## Standalone

Whatever environment (your host or a docker container) you build pyvxl in will require these dependencies:
- g++
- Python 3 development libraries (e.g. `python3-devel`)
- GeoTIFF development libraries (e.g. `libgeotiff-devel`)

If you're on your host, it's recommended you do this in a virtualenv so you don't mess with your

```bash
mkdir build && cd build
cmake -DVXL_DIR=$VXL_SRC_DIR -DPYBIND11_DIR=$PYBIND11_SRC_DIR -G Ninja $PYVXL_SRC_DIR
-DCMAKE_CXX_FLAGS="-fdiagnostics-color=always" -DCMAKE_C_FLAGS="-fdiagnostics-color=always" ..
ninja
ninja install
```

### Using a virtualenv
To avoid messing with your host machine, you can build pyvxl standalone in a virtualenv

```
virtualenv -p python3 env
source env/bin/activate
cmake -DVXL_DIR=${VXL_DIR} -DPYBIND11_DIR=${PYBIND11_DIR} -G Ninja -DCMAKE_CXX_FLAGS="-fdiagnostics-color=always" -DCMAKE_C_FLAGS="-fdiagnostics-color=always" -DPYTHON_SITE=${YOUR_VIRTUALENV_SITE-PACKAGES} ..
ninja
ninja install
```

### Contrib (Optional)

If you want to use any contrib module X, then add PYVXL_CONTRIB_MAKE_X=ON as a CMake option. To build all the contrib packages, use PYVXL_CONTRIB_MAKE_ALL=ON.

Also, uncomment these lines from the main CMakeLists.txt:

```
set(VXL_BUILD_CONTRIB "TRUE" CACHE BOOL "BUILD CONTRIB")
set(VXL_BUILD_CORE_VIDEO "TRUE" CACHE BOOL "BUILD CORE VIDEO")
set(VXL_BUILD_BRL "TRUE" CACHE BOOL "BUILD BRL")
```

### Test your build

```
python3 -c "import vxl"
```

## As part of your project
In your cmake file, just add

```
find_package(PythonLibs 3 REQUIRED)
add_subdirectory(${Your pybind11 dir})
add_subdirectory(${Your pyvxl dir})
```

