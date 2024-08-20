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

# find_package(pybind11)

# find_package(Torch REQUIRED)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

# find_package(PythonLibs REQUIRED)
# include_directories(${PYTHON_INCLUDE_DIRS})
# link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../torchami/install/lib/)
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../torchami/install/include/)



################################
# Tests
################################
# Add test cpp file
#add_executable(blank blank.cpp)
#add_executable(fb fermi_bose.cpp)
#add_executable(helper_functions helper_functions.cpp)
#add_executable(num_test num_test.cpp)
#add_executable(construct_tests construct_tests.cpp)

add_executable(num_test num_test.cpp)
target_link_libraries(num_test gtest gtest_main torchami tamigraph ${PYTHON_LIBRARIES} ${TORCH_LIBRARIES})
add_executable(construct_tests construct_tests.cpp)
target_link_libraries(construct_tests gtest gtest_main torchami tamigraph ${PYTHON_LIBRARIES} ${TORCH_LIBRARIES})

#target_link_libraries(fb torchami tamigraph gtest gtest_main)
#target_link_libraries(helper_functions torchami tamigraph gtest gtest_main)
#target_link_libraries(num_test torchami tamigraph gtest gtest_main ${PYTHON_LIBRARIES} boost_system boost_filesystem)
#target_link_libraries(construct_tests torchami tamigraph gtest gtest_main)

