cmake_minimum_required(VERSION 3.25)
project(Algorithm C)

set(CMAKE_C_STANDARD 23)

add_executable(Algorithm main.c)

########### simd ############
add_subdirectory(simd)
target_link_libraries(Algorithm simd)
