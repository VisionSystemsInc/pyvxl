# pyvxl
Python wrapper for vxl using pybind11

# Requirements

- pybind11 - You need to add the pybind11 repo as a submodule to your proeject (or install it into the OS)


# Cmake

In your cmake file, just add

```
find_package(PythonLibs 3 REQUIRED)
add_subdirectory(${Your pybind11 dir})
add_subdirectory(${Your pyvxl dir})
```
