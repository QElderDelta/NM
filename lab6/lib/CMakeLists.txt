set(CMAKE_CXX_STANDARD 17)

add_library(hyperbolic_dif_eq_solver
        include/hyperbolic_dif_eq_solver.h
        src/hyperbolic_dif_eq_solver.cpp)

target_include_directories(hyperbolic_dif_eq_solver PUBLIC include)
target_link_libraries(hyperbolic_dif_eq_solver PUBLIC MatrixOps)