cmake_minimum_required(VERSION 3.19)
project(lab4)

set(CMAKE_CXX_STANDARD 17)

add_executable(lab4 main.cpp include/cauchyproblemsolvers.h src/cauchyproblemsolvers.cpp ../lab3/include/integral.h
               ../lab3/src/integral.cpp)