# pyvxl
Python wrapper for vxl using pybind11

# Requirements

- pybind11 - You need to add the pybind11 repo as a submodule to your project (or install it into the OS)
- vxl

# Cmake

pyvxl can be build standalone or as a module of an existing project.

## Standalone
```bash
cmake -DVXL_DIR=$VXL_SRC_DIR -DPYBIND11_DIR=$PYBIND11_SRC_DIR $PYVXL_SRC_DIR
```

## As part of your project
In your cmake file, just add

```
find_package(PythonLibs 3 REQUIRED)
add_subdirectory(${Your pybind11 dir})
add_subdirectory(${Your pyvxl dir})
```
