cmake_minimum_required(VERSION 3.20)

# Fetch GoogleTest (header-only usage for most parts, but we’ll link the libs)
include(FetchContent)

# Allow overriding the tag from the CLI if needed: -DGTEST_TAG=v1.13.0
set(GTEST_TAG v1.14.0 CACHE STRING "googletest tag to fetch")

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        ${GTEST_TAG}
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(googletest)

# Avoid overriding parent project’s settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Torch is already found at top-level; require here for clarity
find_package(Torch REQUIRED)

# The libraries tami_graph and tami_base already export their include dirs,
# so we don’t need include_directories() in tests.

enable_testing()
include(GoogleTest)

# Helper function to declare a test target and register it with ctest
function(add_tami_test name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name}
    PRIVATE
      GTest::gtest_main
      tami_graph       # brings in tami_base transitively
      Torch::Torch
  )
  gtest_discover_tests(${name})
endfunction()

# ---- Tests ----
add_tami_test(num_test         ${CMAKE_CURRENT_SOURCE_DIR}/num_test.cpp)
add_tami_test(construct_tests  ${CMAKE_CURRENT_SOURCE_DIR}/construct_tests.cpp)
add_tami_test(helper_tests     ${CMAKE_CURRENT_SOURCE_DIR}/helper_functions.cpp)
add_tami_test(fb               ${CMAKE_CURRENT_SOURCE_DIR}/fermi_bose.cpp)
