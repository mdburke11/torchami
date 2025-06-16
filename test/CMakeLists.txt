cmake_minimum_required(VERSION 3.18.0)


project(googletest-git NONE)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.10.0
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

find_package(Torch REQUIRED)
cmake_print_variables(TORCH_CXX_FLAGS)
cmake_print_variables(TORCH_INCLUDE_DIRS)
include_directories(${TORCH_INCLUDE_DIRS})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

include_directories(${CMAKE_SOURCE_DIR}/src/tami_graph_src)
include_directories(${CMAKE_SOURCE_DIR}/src/tami_base_src)

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

################################
# Tests
################################
# Add test cpp file
add_executable(num_test num_test.cpp)
target_link_libraries(num_test gtest gtest_main torchami tamigraph ${PYTHON_LIBRARIES} ${TORCH_LIBRARIES})
add_executable(construct_tests construct_tests.cpp)
target_link_libraries(construct_tests gtest gtest_main torchami tamigraph ${PYTHON_LIBRARIES} ${TORCH_LIBRARIES})
add_executable(helper_tests helper_functions.cpp)
target_link_libraries(helper_tests gtest gtest_main torchami tamigraph ${PYTHON_LIBRARIES} ${TORCH_LIBRARIES})
add_executable(fb fermi_bose.cpp)
target_link_libraries(fb gtest gtest_main torchami tamigraph ${PYTHON_LIBRARIES} ${TORCH_LIBRARIES})
