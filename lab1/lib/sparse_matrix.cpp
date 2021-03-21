#include "sparse_matrix.h"

std::istream &operator>>(std::istream &os, TSparseMatrix &matrix) {
    size_t n, m;
    double val;
    size_t nnz_count = 0;
    os >> n >> m;
    for(size_t i = 0; i < n; ++i) {
        matrix.rows_counts_.push_back(nnz_count);
        for(size_t j = 0; j < m; ++j) {
            os >> val;
            if(val != 0) {
                matrix.values_.push_back(val);
                matrix.cols_idx_.push_back(j);
                ++nnz_count;
            }
        }
    }
    matrix.rows_counts_.push_back(nnz_count);
    matrix.nnz_count_ = nnz_count;
    return os;
}

//multiplies sparse matrix by a vector
TMatrix operator*(const TSparseMatrix &a, const TMatrix &b) {
    assert(b.GetSize().number_of_cols == 1);
    size_t n = b.GetSize().number_of_rows;
    TMatrix c(b.GetSize());
    c.Clear();
    for(size_t i = 0; i < n; ++i) {
        for(size_t j = a.rows_counts_[i]; j < a.rows_counts_[i + 1]; ++j) {
            c.SetElement(i, 0, c.GetElement(i, 0) + a.values_[j] * b.GetElement(a.cols_idx_[j], 0));
        }
    }
    return c;
}
