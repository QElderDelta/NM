cmake_minimum_required(VERSION 3.21)
project(lab7)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(lib)
add_subdirectory(${CMAKE_SOURCE_DIR}/../MatrixOps cmake-build-debug)

add_executable(lab7 main.cpp)
target_link_libraries(lab7 PUBLIC elliptic_equation_solver)
