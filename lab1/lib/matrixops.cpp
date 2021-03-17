#include "matrixops.h"

TMatrix BackwardSubstitution(const TMatrix &a, const TMatrix &b, size_t column_number) {
    TMatrix x(a.GetSize().number_of_rows, 1);
    double row_sum;
    for(int i = a.GetSize().number_of_rows - 1; i >= 0; --i) {
        row_sum = 0;
        for(size_t j = i + 1; j < a.GetSize().number_of_cols; ++j) {
            row_sum += a.GetElement(i, j) * x.GetElement(j, 0);
        }
        x.SetElement(i, 0, (b.GetElement(i, column_number) - row_sum) / a.GetElement(i, i));
    }
    return x;
}

TMatrix ForwardSubstitution(const TMatrix &a, const TMatrix &b, size_t column_number) {
    TMatrix x(a.GetSize().number_of_rows, 1);
    double row_sum;
    for(size_t i = 0; i < a.GetSize().number_of_rows; ++i) {
        row_sum = 0;
        for(size_t j = 0; j < i; ++j) {
            row_sum += a.GetElement(i, j) * x.GetElement(j, 0);
        }
        x.SetElement(i, 0, b.GetElement(i, column_number) - row_sum);
    }
    return x;
}

double GetDeterminantUsingLU(const TMatrix &l, const TMatrix &u, const TMatrix &p) {
    bool parity_is_odd = FindParityOfPermutation(p);
    double result = 1;
    double product = 1;
    for(int i = 0; i < l.GetSize().number_of_rows; ++i) {
        product *= l.GetElement(i, i);
    }
    result *= product;
    product = 1;
    for(int i = 0; i < l.GetSize().number_of_rows; ++i) {
        product *= u.GetElement(i, i);
    }
    result *= product;
    if(parity_is_odd) {
        result *= -1;
    }
    return result;
}

TMatrix InverseMatrixUsingLU(const TMatrix &l, const TMatrix &u, const TMatrix &p) {
    TMatrix z(l.GetSize());
    TMatrix e(l.GetSize());
    TMatrix x;
    e = p * e;
    for(size_t j = 0; j < l.GetSize().number_of_cols; ++j) {
        x = ForwardSubstitution(l, e, j);
        for(int i = 0; i < l.GetSize().number_of_rows; ++i) {
            z.SetElement(i, j, x.GetElement(i, 0));
        }
    }
    for(size_t j = 0; j < l.GetSize().number_of_cols; ++j) {
        x = BackwardSubstitution(u, z, j);
        for(int i = 0; i < l.GetSize().number_of_rows; ++i) {
            e.SetElement(i, j, x.GetElement(i, 0));
        }
    }
    return e;
}

TMatrix SolveLinearSystemUsingLU(const TMatrix &l, const TMatrix &u, const TMatrix &p, const TMatrix &b) {
    TMatrix z(l.GetSize().number_of_rows, 1);
    TMatrix x(l.GetSize().number_of_rows, 1);
    TMatrix b_permutation;
    b_permutation = p * b;
    z = ForwardSubstitution(l, b_permutation, 0);
    x = BackwardSubstitution(u, z, 0);
    return x;
}

bool FindParityOfPermutation(TMatrix p) {
    size_t permutation_count = 0;
    for(int i = 0; i < p.GetSize().number_of_rows; ++i) {
        if(p.GetElement(i, i) != 1) {
            ++permutation_count;
            p.SetElement(i, i, 1);
            for(int j = i + 1; j < p.GetSize().number_of_rows; ++j) {
                if(p.GetElement(j, i) == 1) {
                    p.SetElement(j, i, 0);
                }
            }
        }
    }
    return permutation_count % 2;
}

TMatrix SweepMethod(const TTridiagonalMatrix &a, const TMatrix &b) {
    size_t n = b.GetSize().number_of_rows;
    std::vector<std::vector<double>> coefs(n - 1);
    TMatrix x(n, 1);
    coefs[0].push_back(-a.GetElement(0, 1) / a.GetElement(0, 0));
    coefs[0].push_back(b.GetElement(0, 0) / a.GetElement(0, 0));
    for(size_t i = 1; i < n - 1; ++i) {
        coefs[i].push_back(-a.GetElement(i, 2) / (a.GetElement(i, 1) + a.GetElement(i, 0) * coefs[i - 1][0]));
        coefs[i].push_back((b.GetElement(i, 0) - a.GetElement(i, 0) * coefs[i - 1][1]) /
                            (a.GetElement(i, 1) + a.GetElement(i, 0) * coefs[i - 1][0]));
    }
    x.SetElement(n - 1, 0, (b.GetElement(n - 1, 0) - a.GetElement(n - 1, 0) * coefs[n - 1 - 1][1]) /
                           (a.GetElement(n - 1, 1) + a.GetElement(n - 1, 0) * coefs[n - 1 - 1][0]));
    for(int i = n - 2; i >= 0; --i) {
        x.SetElement(i, 0, coefs[i][0] * x.GetElement(i + 1, 0) + coefs[i][1]);
    }
    return x;
}
