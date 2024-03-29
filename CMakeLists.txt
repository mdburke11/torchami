cmake_minimum_required(VERSION 3.18...3.22.1)
project( P )
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/common/cmake)

add_subdirectory(${CMAKE_SOURCE_DIR}/src/tami_base_src/ src/tami_base_src/)

add_subdirectory(${CMAKE_SOURCE_DIR}/src/tami_graph_src/ src/tami_graph_src/)

# See if user has pybind11 installed
find_package(pybind11)
if(pybind11_FOUND)
  
  if(MAKE_PYTAMI)

    add_subdirectory(${CMAKE_SOURCE_DIR}/src/pytami_src/ src/pytami_src/)

  else()
    message("PYTAMI is available. But is not being compiled.  To Enable set MAKE_PYTAMI=ON and set PYTHON_LIBRARY_DIR path for installation.")
  ENDIF(MAKE_PYTAMI)

else()
  message("pybind11 was not found.  Python library will not be compiled. Continuing without.")
endif()

if(BUILD_DOC)


message("Adding doc subdirectory")
add_subdirectory("doc")

#########################################################
#Documentation
find_package(Doxygen)

if(DOXYGEN_FOUND)
  # # exclude sqlite code
  # set(DOXYGEN_EXCLUDE_PATTERNS
        # */sqlite3/*
  # )
  message("run 'make DOCS' to generate documentation")
  # # doxygen settings can be set here, prefixed with "DOXYGEN_"
  set(DOXYGEN_EXTRACT_PRIVATE YES)
  set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/doc/doxygen")
  set(DOXYGEN_GENERATE_LATEX YES)
  set(DOXYGEN_GENERATE_XML YES)
  set(DOXYGEN_SORT_MEMBER_DOCS YES)
  set(DOXYGEN_USE_MATHJAX YES)
  # #set(DOXYGEN_USE_MDFILE_AS_MAINPAGE main_page.dox)
  # # this target will only be built if specifically asked to.
  # # run "make LibAMI-docs" to create the doxygen documentation
  # # run "make DOCS" to create the doxygen documentation
  doxygen_add_docs(
    DOCS
    ${PROJECT_SOURCE_DIR}/src
    COMMENT "Generate API-documents for torchami."
  )
endif(DOXYGEN_FOUND)


endif(BUILD_DOC)
