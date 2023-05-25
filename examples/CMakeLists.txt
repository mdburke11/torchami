cmake_minimum_required(VERSION 3.18...3.22.1)
project(example CXX)

add_executable(${PROJECT_NAME}
./src/pt_main.cpp
./src/examples.hpp
)


# target_link_libraries(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib/libami.so m)
link_directories(${Boost_LIBRARY_DIRS})
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libami/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libami/install/include)

# Changed to compile with c++14 compiler - Maybe can stay as c++11 since these are compiled separately
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -stdlib=libc++")
else()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -DUSE_MPI -lstdc++fs ") 
endif()


target_link_libraries(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/../libami/install/lib/libami.dylib
 m  )
