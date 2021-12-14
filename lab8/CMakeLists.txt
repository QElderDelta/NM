cmake_minimum_required(VERSION 3.21)
project(lab8)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(lib)
add_subdirectory(${CMAKE_SOURCE_DIR}/../MatrixOps cmake-build-debug)

add_executable(lab8 main.cpp)
target_link_libraries(lab8 PUBLIC 2d_equation_solver)
