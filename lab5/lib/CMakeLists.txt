add_library(parabolic_dif_eq_solver
        include/parabolic_dif_eq_solver.h
        src/parabolic_dif_eq_solver.cpp)

target_include_directories(parabolic_dif_eq_solver PUBLIC include)
target_link_libraries(parabolic_dif_eq_solver PUBLIC MatrixOps)