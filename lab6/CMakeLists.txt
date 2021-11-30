cmake_minimum_required(VERSION 3.20)
project(lab6)

add_subdirectory(lib)
add_subdirectory(${CMAKE_SOURCE_DIR}/../MatrixOps cmake-build-debug)
set(CMAKE_CXX_STANDARD 17)

add_executable(lab6 main.cpp)
target_link_libraries(lab6 PUBLIC hyperbolic_dif_eq_solver)