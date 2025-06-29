cmake_minimum_required(VERSION 3.18...3.22.1)
project( tamigraph )
add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)

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

set(LIB_GRAPH_FILES
  ${CMAKE_SOURCE_DIR}/src/tami_graph_src/tamigraph_helpers.cpp
  ${CMAKE_SOURCE_DIR}/src/tami_graph_src/tamigraph_labelling.cpp
  ${CMAKE_SOURCE_DIR}/src/tami_graph_src/tamigraph_interface.cpp
  ${CMAKE_SOURCE_DIR}/src/tami_graph_src/tamigraph_constructs.cpp
)

ADD_LIBRARY(tamigraph SHARED ${LIB_GRAPH_FILES})
include_directories(${CMAKE_SOURCE_DIR}/src/tami_graph_src)
include_directories(${CMAKE_SOURCE_DIR}/src/tami_base_src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../tami_base_src)

#link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../install/lib)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../install/include/tami_base_src)

#message(STATUS "libtorchami directory: ${CMAKE_CURRENT_SOURCE_DIR}/../../install/lib")

#target_link_libraries(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/../../install/lib/libtorchami.so
# m )


target_link_libraries(${PROJECT_NAME} m ${TORCH_LIBRARIES})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} m ${PYTHON_LIBRARIES})

INSTALL (TARGETS tamigraph
ARCHIVE DESTINATION "lib/"
LIBRARY DESTINATION "lib/"
COMPONENT library)

# INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/src/."
#         DESTINATION include
#         FILES_MATCHING PATTERN "*.hpp"
# )

# Preserve tami_base_src and tami_graph_src folders in include path

install(DIRECTORY "${CMAKE_SOURCE_DIR}/src/tami_graph_src"
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp")

