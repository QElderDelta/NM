#include "tridiagonal_matrix.h"

std::istream &operator>>(std::istream &is, TTridiagonalMatrix &m) {
    is >> m.number_of_rows_;
    m.matrix_.resize(m.number_of_rows_);
    m.matrix_[0].resize(2);
    is >> m.matrix_[0][0] >> m.matrix_[0][1];
    for(size_t i = 1; i < m.number_of_rows_ - 1; ++i) {
        m.matrix_[i].resize(3);
        for(size_t j = 0; j < 3; ++j) {
            is >> m.matrix_[i][j];
        }
    }
    m.matrix_[m.number_of_rows_ - 1].resize(2);
    is >> m.matrix_[m.number_of_rows_ - 1][0] >> m.matrix_[m.number_of_rows_ - 1][1];
    return is;
}

double &TTridiagonalMatrix::GetElement(size_t i, size_t j) {
    if(i > number_of_rows_ || ((i == 0 || i == number_of_rows_ - 1) && j > 1) || j > 2 ) {
        throw std::out_of_range("");
    }
    return matrix_[i][j];
}

const double &TTridiagonalMatrix::GetElement(size_t i, size_t j) const {
    if(i > number_of_rows_ || ((i == 0 || i == number_of_rows_ - 1) && j > 1) || j > 2 ) {
        throw std::out_of_range("");
    }
    return matrix_[i][j];
}

void TTridiagonalMatrix::SetElement(size_t i, size_t j, double value) {
    if(i > number_of_rows_ || ((i == 0 || i == number_of_rows_ - 1) && j > 1) || j > 2 ) {
        throw std::out_of_range("");
    }
    matrix_[i][j] = value;
}
