#include <fstream>
#include <iostream>
#include "lib/matrix.h"

void Task1_1() {
    TMatrix a, b;
    std::ifstream is("../task1_1.txt");
    std::ofstream os("../task1_1_result.txt");
    is >> a >> b;
    auto [l, u, p] = a.LUDecomposition();
    os << "Matrix L:" << '\n' << l << '\n';
    os << "Matrix U:" << '\n' << u << '\n';
    os << "P * L * U:" << '\n' << p * l * u << '\n';
    os << "Ax = b solution:" << '\n' << SolveLinearSystemUsingLU(l, u, p, b) << '\n';
    os << "Determinant of matrix A: " << GetDeterminantUsingLU(l, u, p) << '\n';
    os << "Inverse of matrix A:" << '\n' << InverseMatrixUsingLU(l, u, p) << '\n';
    os << "Product of matrix A and its' inverse matrix" << '\n'<< a * InverseMatrixUsingLU(l, u, p) << '\n';
}

int main() {
    Task1_1();
    return 0;
}
