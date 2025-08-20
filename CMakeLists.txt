cmake_minimum_required(VERSION 3.20)
project(torchami LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(GNUInstallDirs)

option(MAKE_PYTAMI "Build the pybind11 Python extension" OFF)
option(BUILD_DOC "Build documentation" OFF)

# Always need a Python interpreter (used to query Torch info, etc.)
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# ---------------------------------------------------------------------------
# Torch
# ---------------------------------------------------------------------------
find_package(Torch QUIET)
if(NOT Torch_FOUND)
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import torch,sys;print(torch.utils.cmake_prefix_path)"
    OUTPUT_VARIABLE _torch_prefix
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(_torch_prefix)
    list(APPEND CMAKE_PREFIX_PATH "${_torch_prefix}")
    find_package(Torch REQUIRED)
  else()
    message(FATAL_ERROR "PyTorch not found and torch.utils.cmake_prefix_path not available.")
  endif()
endif()

# Match ABI flag to the installed PyTorch
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c "import torch;print(int(torch._C._GLIBCXX_USE_CXX11_ABI))"
  OUTPUT_VARIABLE _torch_abi
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(_torch_abi MATCHES "^[01]$")
  add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=${_torch_abi})
else()
  add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)
endif()

# Ensure an imported target exists even on legacy Torch exports
if(NOT TARGET Torch::Torch)
  message(STATUS "Torch::Torch target not provided by this Torch package; creating a compatibility target")
  add_library(Torch::Torch INTERFACE IMPORTED)
  if(DEFINED TORCH_INCLUDE_DIRS)
    target_include_directories(Torch::Torch INTERFACE "${TORCH_INCLUDE_DIRS}")
  endif()
  if(DEFINED TORCH_CXX_FLAGS AND TORCH_CXX_FLAGS)
    separate_arguments(_TORCH_CXX_FLAGS NATIVE_COMMAND "${TORCH_CXX_FLAGS}")
    target_compile_options(Torch::Torch INTERFACE ${_TORCH_CXX_FLAGS})
  endif()
  if(DEFINED TORCH_LIBRARIES)
    target_link_libraries(Torch::Torch INTERFACE ${TORCH_LIBRARIES})
  endif()
endif()

# ---------------------------------------------------------------------------
# pybind11 / Python dev (ONLY if building the Python module)
# ---------------------------------------------------------------------------
if(MAKE_PYTAMI)
  # Need dev headers + module support for pybind11 new tools
  find_package(Python3 COMPONENTS Interpreter Development Development.Module REQUIRED)

  # Ask pybind11 to use CMake's FindPython if it can
  set(PYBIND11_FINDPYTHON ON)

  # Prefer CONFIG packages (pip installs provide this)
  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    # Autodetect pybind11 config dir from Python
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
      OUTPUT_VARIABLE _pybind11_dir
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    if(_pybind11_dir AND EXISTS "${_pybind11_dir}/pybind11Config.cmake")
      list(APPEND CMAKE_PREFIX_PATH "${_pybind11_dir}")
      find_package(pybind11 CONFIG REQUIRED)
    else()
      message(FATAL_ERROR "MAKE_PYTAMI=ON but pybind11 not found; install it (pip) or set pybind11_DIR.")
    endif()
  endif()
endif()

# ---------------------------------------------------------------------------
# Boost
# ---------------------------------------------------------------------------
find_package(Boost 1.70 REQUIRED)

# ---------------------------------------------------------------------------
# Subdirectories
# ---------------------------------------------------------------------------
add_subdirectory(src/tami_base_src)
add_subdirectory(src/tami_graph_src)

if(MAKE_PYTAMI)
  add_subdirectory(src/pytami_src)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(CMAKE_BUILD_TYPE MATCHES "Release|RelWithDebInfo|Debug" OR DEFINED TEST)
  if(EXISTS "${CMAKE_SOURCE_DIR}/test/CMakeLists.txt")
    message(STATUS "Testing enabled")
    enable_testing()
    add_subdirectory(test)
  endif()
endif()

# ---------------------------------------------------------------------------
# Docs
# ---------------------------------------------------------------------------
if(BUILD_DOC)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_EXTRACT_PRIVATE YES)
    set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/doc/doxygen")
    set(DOXYGEN_GENERATE_LATEX YES)
    set(DOXYGEN_GENERATE_XML YES)
    set(DOXYGEN_SORT_MEMBER_DOCS YES)
    set(DOXYGEN_USE_MATHJAX YES)
    doxygen_add_docs(
      DOCS
      ${PROJECT_SOURCE_DIR}/src
      COMMENT "Generate API documents for torchami."
    )
  endif()
endif()

# ---- (Optional) debug prints ----
message(STATUS "Python3_VERSION            = ${Python3_VERSION}")
message(STATUS "Python3_EXECUTABLE         = ${Python3_EXECUTABLE}")
message(STATUS "Torch_DIR                  = ${Torch_DIR}")
message(STATUS "TORCH_INCLUDE_DIRS         = ${TORCH_INCLUDE_DIRS}")
message(STATUS "TORCH_LIBRARIES            = ${TORCH_LIBRARIES}")
message(STATUS "TORCH_CXX_FLAGS            = ${TORCH_CXX_FLAGS}")
