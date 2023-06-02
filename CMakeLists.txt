cmake_minimum_required(VERSION 3.18...3.22.1)
project(torchami CXX)


find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

link_directories(${Boost_LIBRARY_DIRS})

# Changed to compile with c++14 compiler - Maybe can stay as c++11 since these are compiled separately
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -stdlib=libc++")
else()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -lstdc++fs ") 
endif()

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -lstdc++fs")# -msse4.2")# -ffast-math -funsafe-math-optimizations -msse4.2")# -ffast-math")

set(LIB_FILES
  ./src/tami_ft.cpp
  ./src/tami_construct.cpp
  ./src/tami_evaluate.cpp
  ./src/tami_terms.cpp
)

ADD_LIBRARY(torchami SHARED ${LIB_FILES})
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

#install(TARGETS dmft RUNTIME DESTINATION "bin/")

