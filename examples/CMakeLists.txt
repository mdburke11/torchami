cmake_minimum_required(VERSION 3.18...3.22.1)
project(example CXX)

add_executable(${PROJECT_NAME}
./src/pt_main.cpp
./src/examples.hpp
./src/examples.cpp
)

find_package(Boost REQUIRED)


link_directories(${Boost_LIBRARY_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../torchami/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../torchami/install/include)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -stdlib=libc++")
else()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -lstdc++fs ") 
endif()

if (MACOSX)
    target_link_libraries(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/../../torchami/install/lib/libtorchami.dylib
    m  )
else()
    target_link_libraries(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/../../torchami/install/lib/libtorchami.so
    m  )
endif()
