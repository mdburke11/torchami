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

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})


################################
# Tests
################################
# Add test cpp file
#add_executable(blank blank.cpp)
#add_executable(fb fermi_bose.cpp)
#add_executable(helper_functions helper_functions.cpp)
add_executable(num_test num_test.cpp)
#add_executable(construct_tests construct_tests.cpp)

#target_link_libraries(fb torchami tamigraph gtest gtest_main)
#target_link_libraries(helper_functions torchami tamigraph gtest gtest_main)
target_link_libraries(num_test torchami tamigraph gtest gtest_main ${PYTHON_LIBRARIES})
#target_link_libraries(construct_tests torchami tamigraph gtest gtest_main)

