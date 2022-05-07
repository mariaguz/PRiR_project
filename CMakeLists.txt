cmake_minimum_required(VERSION 3.15)
project(lab_2)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP REQUIRED)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(lab_2 main.cpp)
