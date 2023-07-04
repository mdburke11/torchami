cmake_minimum_required(VERSION 3.18...3.22.1)
project( tamigraph )

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

link_directories(${Boost_LIBRARY_DIRS})

set(LIB_BASE_FILES
  ${CMAKE_SOURCE_DIR}/src/tami_graph_src/main.cpp
)

INSTALL (TARGETS torchami 
ARCHIVE DESTINATION "lib/"
LIBRARY DESTINATION "lib/"
COMPONENT library)

INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/src/."
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
)