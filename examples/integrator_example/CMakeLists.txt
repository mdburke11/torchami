cmake_minimum_required(VERSION 3.16...3.22.1)
project(integrator_example CXX)

add_executable(${PROJECT_NAME}
./src/integration.hpp
./src/graph_util.cpp
./src/integrator_eg.cpp
./src/main.cpp
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -stdlib=libc++")
else()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -lstdc++fs ") 
endif()

find_package(Boost REQUIRED)
link_directories(${Boost_LIBRARY_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
target_link_libraries(${PROJECT_NAME} ${TORCH_LIBRARIES})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} m ${PYTHON_LIBRARIES})

if(LIBTORCHAMI_ROOT)
link_directories(${LIBTORCHAMI_ROOT}/lib/tami_base_src)
include_directories(${LIBTORCHAMI_ROOT}/include/tami_base_src)
link_directories(${LIBTORCHAMI_ROOT}/lib/tami_graph_src)
include_directories(${LIBTORCHAMI_ROOT}/include/tami_graph_src)

IF(MACOSX)
set(CMAKE_MACOSX_RPATH=OFF)
target_link_libraries(${PROJECT_NAME}  ${LIBTORCHAMI_ROOT}/lib/libtorchami.dylib  ${LIBTORCHAMI_ROOT}/lib/libtamigraph.dylib
 m )
 else()
	 target_link_libraries(${PROJECT_NAME}  ${LIBTORCHAMI_ROOT}/lib/libtorchami.so ${LIBTORCHAMI_ROOT}/lib/libtamigraph.so
 m )
 ENDIF(MACOSX)
 
else()
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../torchami/install/lib/tami_base_src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../torchami/install/include/tami_base_src)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../torchami/install/lib/tami_graph_src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../torchami/install/include/tami_graph_src)

IF(MACOSX)
set(CMAKE_MACOSX_RPATH=OFF)
target_link_libraries(${PROJECT_NAME}  ${CMAKE_CURRENT_SOURCE_DIR}/../../../torchami/install/lib/libtorchami.dylib ${CMAKE_CURRENT_SOURCE_DIR}/../../../torchami/install/lib/libtamigraph.dylib
 m )
 else()
target_link_libraries(${PROJECT_NAME}  ${CMAKE_CURRENT_SOURCE_DIR}/../../../torchami/install/lib/libtorchami.so ${CMAKE_CURRENT_SOURCE_DIR}/../../../torchami/install/lib/libtamigraph.so
 m )
 ENDIF(MACOSX)
ENDIF(LIBTORCHAMI_ROOT)
