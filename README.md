# pyvxl
Python wrapper for vxl using pybind11

[![CircleCI](https://circleci.com/gh/VisionSystemsInc/pyvxl.svg?style=svg)](https://circleci.com/gh/VisionSystemsInc/pyvxl)

# Requirements

- pybind11 - You need to add the pybind11 repo as a submodule to your project (or install it into the OS)
- vxl

# Cmake

pyvxl can be build standalone or as a module of an existing project.

## Standalone

Whatever environment (your host or a docker container) you build pyvxl in will require these dependencies
- g++
- Python 3 development libraries (e.g. `python3-devel`)
- GeoTIFF development libraries (e.g. `libgeotiff-devel`)
- CMake 3.10.2 or newer
- Ninja (optional)

If you're not using ninja, replace "ninja" in the commands below with "make".

### Using a virtualenv (recommended)

To avoid messing with your host machine, you can build pyvxl in a virtualenv

```
virtualenv -p python3 env
source env/bin/activate
mkdir build && cd build
cmake -DVXL_DIR=${VXL_DIR} -DPYBIND11_DIR=${PYBIND11_DIR} -G Ninja -DPYTHON_SITE=${YOUR_VIRTUALENV_SITE-PACKAGES} ..
ninja
ninja install
```

### Not using virtualenv

If you're not going to use a virtualenv because you're in a docker, or your host machine (not recommended), then you can drop the PYTHON_SITE variable

```bash
mkdir build && cd build
cmake -DVXL_DIR=${VXL_DIR} -DPYBIND11_DIR=${PYBIND11_DIR} -G Ninja ..
ninja
ninja install
```

### Contrib (Optional)

If you want to use any contrib module X, then add PYVXL_CONTRIB_MAKE_X=ON as a CMake option. To build all the contrib packages, use PYVXL_CONTRIB_MAKE_ALL=ON.

Also, these VXL build options should be set:

```
-DVXL_BUILD_CONTRIB="TRUE" -DVXL_BUILD_CORE_VIDEO="TRUE" -DVXL_BUILD_BRL="TRUE" 
```

### Test your build

```
python3 -c "import vxl"
```

## As part of your project
In your CMake file, just add

```
find_package(PythonLibs 3 REQUIRED)
add_subdirectory(${Your pybind11 dir})
add_subdirectory(${Your pyvxl dir})
```

