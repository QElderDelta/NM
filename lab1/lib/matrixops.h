#pragma once

#include "matrix.h"
#include "tridiagonal_matrix.h"

class TMatrix;

TMatrix BackwardSubstitution(const TMatrix& a, const TMatrix& b, size_t column_number);

bool FindParityOfPermutation(TMatrix p);

TMatrix ForwardSubstitution(const TMatrix& a, const TMatrix& b, size_t column_number);

double GetDeterminantUsingLU(const TMatrix& l, const TMatrix& u, const TMatrix& p);

TMatrix InverseMatrixUsingLU(const TMatrix& l, const TMatrix& u, const TMatrix& p);

TMatrix SeidelMethod(const TMatrix& m, const TMatrix& b, double eps, std::ostream& log_stream = std::cout);

TMatrix SimpleIterationsMethod(const TMatrix& m, const TMatrix& b, double eps, std::ostream& log_stream = std::cout);

TMatrix SolveLinearSystemUsingLU(const TMatrix& l, const TMatrix& u, const TMatrix& p, const TMatrix& b);

TMatrix SweepMethod(const TTridiagonalMatrix& a, const TMatrix& b);