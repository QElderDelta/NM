set(CMAKE_CXX_STANDARD 17)

add_library(2d_equation_solver include/2d_eq_solver.h src/2d_eq_solver.cpp)

target_include_directories(2d_equation_solver PUBLIC include)
target_link_libraries(2d_equation_solver PUBLIC MatrixOps)