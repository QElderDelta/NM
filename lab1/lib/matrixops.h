#pragma once

#include <complex>
#include <variant>

#include "matrix.h"
#include "tridiagonal_matrix.h"

class TMatrix;

using Eigenvalue = std::variant<double, std::complex<double>>;

TMatrix BackwardSubstitution(const TMatrix& a, const TMatrix& b, size_t column_number);

bool FindParityOfPermutation(TMatrix p);

TMatrix ForwardSubstitution(const TMatrix& a, const TMatrix& b, size_t column_number);

double GetDeterminantUsingLU(const TMatrix& l, const TMatrix& u, const TMatrix& p);

std::vector<std::complex<double>> GetEigenvaluesUsingQR(const TMatrix& a, double eps,
                                                        std::ostream& log_stream = std::cout);

double GetSumOfSquaredNonDiagonalElements(const TMatrix& a);

TMatrix InverseMatrixUsingLU(const TMatrix& l, const TMatrix& u, const TMatrix& p);

std::pair<TMatrix, TMatrix> JacobiRotationMethod(const TMatrix& a, double eps, std::ostream& log_stream = std::cout);

TMatrix SeidelMethod(const TMatrix& m, const TMatrix& b, double eps, std::ostream& log_stream = std::cout);

TMatrix SimpleIterationsMethod(const TMatrix& m, const TMatrix& b, double eps, std::ostream& log_stream = std::cout);

TMatrix SolveLinearSystemUsingLU(const TMatrix& l, const TMatrix& u, const TMatrix& p, const TMatrix& b);

std::vector<std::complex<double>> SolveQudraticEquationForQR(const TMatrix& a, size_t column);

TMatrix SweepMethod(const TTridiagonalMatrix& a, const TMatrix& b);