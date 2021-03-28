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
    os << "Ax = b solution:" << '\n';
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
    os << "Ax = b solution:" << '\n';
    os << SeidelMethod(a, b, eps, os) << '\n';
    os << "Exact solution:" << '\n';
    auto [l, u, p] = a.LUDecomposition();
    os << SolveLinearSystemUsingLU(l, u, p, b) << '\n';
}

void Task1_4() {
    TMatrix a;
    double eps;
    std::ifstream is("../task1_4.txt");
    std::ofstream os("../task1_4_result.txt");
    is >> a >> eps;
    auto [aa, u] = JacobiRotationMethod(a, eps, os);
    os << "Eigenvalues:" << '\n';
    for(size_t i = 0; i < a.GetSize().number_of_rows; ++i) {
        os << "a_" << i + 1 << " = " << aa.GetElement(i, i) << '\n';
    }
    os << "Eigenvectors:" << '\n';
    for(size_t i = 0; i < a.GetSize().number_of_rows; ++i) {
        os << "x_" << i + 1 << " = " << "(";
        for(size_t j = 0; j < a.GetSize().number_of_rows; ++j) {
            os << u.GetElement(i, j);
            if(j != a.GetSize().number_of_rows - 1) {
                os << ", ";
            }
        }
        os << ")" << '\n';
    }
    TMatrix x(a.GetSize().number_of_rows, 1);
    for(size_t i = 0; i < a.GetSize().number_of_rows; ++i) {
        for(size_t j = 0; j < a.GetSize().number_of_rows; ++j) {
            x.SetElement(j, 0, u.GetElement(j, i));
        }
        os << "A * x_" << i + 1 << ":" << '\n' << a * x << '\n';
        os << "a_" << i + 1 << " * x_" << i + 1 << ":" << '\n' << x * aa.GetElement(i, i) << '\n';
        x.Clear();
    }
}

void Task1_5() {
    TMatrix a;
    double eps;
    std::ifstream is("../task1_4.txt");
    std::ofstream os("../task1_5_result.txt");
    is >> a >> eps;
    auto [q, r] = a.QRDecomposition();
    std::cout << q * r << '\n';
}

int main() {
    Task1_1();
    Task1_2();
    Task1_3_1();
    Task1_3_2();
    Task1_4();
    Task1_5();
    return 0;
}
