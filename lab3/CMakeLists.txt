cmake_minimum_required(VERSION 3.19)
project(lab3)

set(CMAKE_CXX_STANDARD 17)

add_executable(lab3 main.cpp include/interpolators.h src/interpolators.cpp ../lab1/lib/matrix.h
               ../lab1/lib/matrix.cpp ../lab1/lib/matrixops.h ../lab1/lib/matrixops.cpp ../lab1/lib/tridiagonal_matrix.h
               ../lab1/lib/tridiagonal_matrix.cpp)