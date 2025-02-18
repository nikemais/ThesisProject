cmake_minimum_required(VERSION 3.10)
project(CPP_SCRATCHPAD)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)
include_directories(${CMAKE_SOURCE_DIR}/include/eigen3)

add_executable(scratchpad
    src/main.cpp 
    src/shadow_calculator.cpp 
    src/utils.cpp 
    src/geometry.cpp
)
