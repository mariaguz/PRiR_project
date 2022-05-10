cmake_minimum_required(VERSION 3.15)
project(projekt)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP REQUIRED)

find_package(MPI REQUIRED)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(projekt main.cpp)

target_link_libraries(projekt ${MPI_C_LIBRARIES})
