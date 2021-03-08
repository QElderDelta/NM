cmake_minimum_required(VERSION 3.17)
project(lab1)

set(CMAKE_CXX_STANDARD 20)

add_executable(lab1 main.cpp lib/matrix.h lib/matrix.cpp lib/matrixops.h lib/matrixops.cpp)