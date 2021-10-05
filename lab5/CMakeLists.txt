cmake_minimum_required(VERSION 3.20)
project(lab5)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(${CMAKE_SOURCE_DIR}/../MatrixOps cmake-build-debug)
add_subdirectory(lib)

add_executable(lab5 main.cpp)
target_link_libraries(lab5 PUBLIC parabolic_dif_eq_solver)