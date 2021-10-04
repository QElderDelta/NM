#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include "matrix.h"

class TSparseMatrix {
public:
    TSparseMatrix() = default;

    friend std::istream& operator>>(std::istream& os, TSparseMatrix& m);

    friend TMatrix operator*(const TSparseMatrix& a, const TMatrix& b);
private:
    size_t nnz_count_{};
    std::vector<double> values_;
    std::vector<size_t> cols_idx_;
    std::vector<size_t> rows_counts_;
};