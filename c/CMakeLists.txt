cmake_minimum_required(VERSION 2.8)
project(SHScont)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11 ")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp ")

add_subdirectory(/home/anthony/libraries/dlib dlib_build)
add_subdirectory(potentials)
include_directories(potentials)
add_subdirectory(opt)
include_directories(opt)
add_subdirectory(exe)

file(COPY input DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/exe)

