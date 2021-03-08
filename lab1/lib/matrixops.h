#pragma once

#include "matrix.h"

class TMatrix;

TMatrix BackwardSubstitution(const TMatrix& a, const TMatrix& b, size_t column_number);

bool FindParityOfPermutation(TMatrix p);

TMatrix ForwardSubstitution(const TMatrix& a, const TMatrix& b, size_t column_number);

double GetDeterminantUsingLU(const TMatrix& l, const TMatrix& u, const TMatrix& p);

TMatrix InverseMatrixUsingLU(const TMatrix& l, const TMatrix& u, const TMatrix& p);

TMatrix SolveLinearSystemUsingLU(const TMatrix& l, const TMatrix& u, const TMatrix& p, const TMatrix& b);