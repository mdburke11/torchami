cmake_minimum_required(VERSION 3.18...3.22.1)
project( pytami )

set(CMAKE_CXX_RELEASE "-O3") # not sure if needed. It was in the example
set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ version selection")

include_directories("${CMAKE_SOURCE_DIR}/src")
include_directories("${CMAKE_SOURCE_DIR}/pytami_src")

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

find_library(TORCH_PYTHON_LIBRARY torch_python PATHS "${TORCH_INSTALL_PREFIX}/lib")
message(STATUS "TORCH_PYTHON_LIBRARY: ${TORCH_PYTHON_LIBRARY}")

file (GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/tami_base_src/*.cpp")
file (GLOB HEADER_FILES "${CMAKE_SOURCE_DIR}/src/tami_base_src/*.hpp")
file (GLOB PYTHON_FILES "${CMAKE_SOURCE_DIR}/pytami_src/*.cpp" "${CMAKE_SOURCE_DIR}/pytami_src/*.hpp")

pybind11_add_module(pytami
    ${SOURCE_FILES}
    ${HEADER_FILES}
    ${PYTHON_FILES}
)

#target_link_libraries(pytami PUBLIC)
target_link_libraries(pytami PUBLIC "${TORCH_LIBRARIES}" ${TORCH_PYTHON_LIBRARY})

install(TARGETS pytami
    COMPONENT python
    LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
)