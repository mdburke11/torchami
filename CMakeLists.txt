cmake_minimum_required(VERSION 3.18...3.22.1)
project( P )

add_subdirectory(${CMAKE_SOURCE_DIR}/src/tami_base_src/)

# See if user has pybind11 installed
find_package(pybind11)
if(pybind11_FOUND)
  
  if(MAKE_PYTAMI)

    add_subdirectory(${CMAKE_SOURCE_DIR}/src/pytami_src/)

  else()
    message("PYTAMI is available. But is not being compiled.  To Enable set MAKE_PYTAMI=ON and set PYTHON_LIBRARY_DIR path for installation.")
  ENDIF(MAKE_PYTAMI)

else()
  message("pybind11 was not found.  Python library will not be compiled. Continuing without.")
endif()
