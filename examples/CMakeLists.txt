cmake_minimum_required(VERSION 3.18...3.22.1)
project(example CXX)

add_executable(${PROJECT_NAME}
./src/pt_main.cpp
./src/examples.hpp
./src/examples.cpp
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
link_directories(${LIBTORCHAMI_ROOT}/lib)
include_directories(${LIBTORCHAMI_ROOT}/include)


IF(MACOSX)
set(CMAKE_MACOSX_RPATH=OFF)
target_link_libraries(${PROJECT_NAME}  ${LIBTORCHAMI_ROOT}/lib/libtorchami.dylib
 m )
 else()
target_link_libraries(${PROJECT_NAME}  ${LIBAMI_ROOT}/lib/libtorchami.so
 m )./else()
 ENDIF(MACOSX)
 
else()
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../torchami/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../torchami/install/include)

IF(MACOSX)
set(CMAKE_MACOSX_RPATH=OFF)
target_link_libraries(${PROJECT_NAME}  ${CMAKE_CURRENT_SOURCE_DIR}/../../torchami/install/lib/libtorchami.dylib
 m )
 else()
target_link_libraries(${PROJECT_NAME}  ${CMAKE_CURRENT_SOURCE_DIR}/../../torchami/install/lib/libtorchami.so
 m )
 ENDIF(MACOSX)
ENDIF(LIBTORCHAMI_ROOT)
