#include <fstream>
#include <iostream>

#include "lib/matrix.h"
#include "lib/sparse_matrix.h"
#include "lib/tridiagonal_matrix.h"

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

void Task1_2() {
    TTridiagonalMatrix a;
    TMatrix b;
    std::ifstream is("../task1_2.txt");
    std::ofstream os("../task1_2_result.txt");
    is >> a >> b;
    os << "Ax = b solution:";
    os << '\n' << SweepMethod(a, b) << '\n';
}

void Task1_3_1() {
    TMatrix a, b;
    double eps;
    std::ifstream is("../task1_3.txt");
    std::ofstream os("../task1_3_1_result.txt");
    is >> a >> b >> eps;
    os << SimpleIterationsMethod(a, b, eps, os) << '\n';
    os << "Exact solution:" << '\n';
    auto [l, u, p] = a.LUDecomposition();
    os << SolveLinearSystemUsingLU(l, u, p, b) << '\n';
}

void Task1_3_2() {
    TMatrix a, b;
    double eps;
    std::ifstream is("../task1_3.txt");
    std::ofstream os("../task1_3_2_result.txt");
    is >> a >> b >> eps;
    os << SeidelMethod(a, b, eps, os) << '\n';
    os << "Exact solution:" << '\n';
    auto [l, u, p] = a.LUDecomposition();
    os << SolveLinearSystemUsingLU(l, u, p, b) << '\n';
}

int main() {
    Task1_1();
    Task1_2();
    Task1_3_1();
    Task1_3_2();
    return 0;
}
