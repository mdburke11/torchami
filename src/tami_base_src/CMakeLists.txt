# Library: tami_base
add_library(tami_base)

# Sources (list explicitly or glob within this dir)
target_sources(tami_base PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/tami_ft.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tami_construct.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tami_evaluate.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tami_terms.cpp
)

# Public headers live here
target_include_directories(tami_base
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link deps
target_link_libraries(tami_base
  PRIVATE
    Torch::Torch
    Boost::boost
)


include(CheckCXXSourceCompiles)

# Build the list of includes that tami_base will see
get_target_property(_TAMI_BASE_INCS tami_base INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(_TAMI_BASE_PRIV  tami_base INCLUDE_DIRECTORIES)
set(_ALL_INCS ${_TAMI_BASE_INCS} ${_TAMI_BASE_PRIV} ${Python3_INCLUDE_DIRS})

# Temporarily set the required includes for the check
set(CMAKE_REQUIRED_INCLUDES "${_ALL_INCS}")
check_cxx_source_compiles("#include <Python.h>\nint main(){return 0;}" HAVE_PYTHON_H)
message(STATUS "tami_base: HAVE_PYTHON_H = ${HAVE_PYTHON_H}")
unset(CMAKE_REQUIRED_INCLUDES)


# Install library + headers (preserve folder name)
install(TARGETS tami_base
  EXPORT torchamiTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tami_base_src
  FILES_MATCHING PATTERN "*.hpp"
)
