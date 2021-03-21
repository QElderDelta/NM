cmake_minimum_required(VERSION 3.17)
project(lab1)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fsanitize=address")

add_executable(lab1 main.cpp lib/matrix.h lib/matrix.cpp lib/matrixops.h lib/matrixops.cpp lib/tridiagonal_matrix.h lib/tridiagonal_matrix.cpp lib/sparse_matrix.h lib/sparse_matrix.cpp)