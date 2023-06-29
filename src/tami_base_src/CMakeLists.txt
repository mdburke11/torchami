cmake_minimum_required(VERSION 3.18...3.22.1)
project( torchami )

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -stdlib=libc++")
else()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -lstdc++fs ") 
endif()

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

link_directories(${Boost_LIBRARY_DIRS})

set(LIB_BASE_FILES
  ${CMAKE_SOURCE_DIR}/src/tami_base_src/tami_ft.cpp
  ${CMAKE_SOURCE_DIR}/src/tami_base_src/tami_construct.cpp
  ${CMAKE_SOURCE_DIR}/src/tami_base_src/tami_evaluate.cpp
  ${CMAKE_SOURCE_DIR}/src/tami_base_src/tami_terms.cpp
  ${CMAKE_SOURCE_DIR}/src/tami_base_src/tami_python.cpp
)

ADD_LIBRARY(torchami SHARED ${LIB_BASE_FILES})
target_link_libraries(${PROJECT_NAME} m ${TORCH_LIBRARIES})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} m ${PYTHON_LIBRARIES})

INSTALL (TARGETS torchami 
ARCHIVE DESTINATION "lib/"
LIBRARY DESTINATION "lib/"
COMPONENT library)

INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/src/."
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
)