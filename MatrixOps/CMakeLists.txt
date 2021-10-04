add_library(MatrixOps
        include/matrix.h
        include/matrixops.h
        include/tridiagonal_matrix.h
        include/sparse_matrix.h
        src/matrix.cpp
        src/matrixops.cpp
        src/tridiagonal_matrix.cpp
        src/sparse_matrix.cpp)

target_include_directories(MatrixOps PUBLIC include)