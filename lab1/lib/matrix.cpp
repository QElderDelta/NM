#include "matrix.h"

bool TMatrixSize::operator==(const TMatrixSize &other) const {
    return this->number_of_rows == other.number_of_rows && this->number_of_cols == other.number_of_cols;
}

bool TMatrixSize::operator!=(const TMatrixSize &other) const {
    return !(*this == other);
}

bool AreCompatible(const TMatrixSize& first, const TMatrixSize& second) {
    return first.number_of_cols == second.number_of_rows;
}

TMatrix::TMatrix(size_t number_of_rows, size_t number_of_cols) : size_(TMatrixSize{number_of_rows, number_of_cols}) {
    InitializeMatrix();
}

TMatrix::TMatrix(const TMatrixSize &size) : size_(size) {
    InitializeMatrix();
}

std::istream &operator>>(std::istream &is, TMatrix &m) {
    is >> m.size_.number_of_rows >> m.size_.number_of_cols;
    m.InitializeMatrix();
    for(size_t i = 0; i < m.size_.number_of_rows; ++i) {
        for(size_t j = 0; j < m.size_.number_of_cols; ++j) {
            is >> m.matrix_[i][j];
        }
    }
    return is;
}

std::ostream &operator<<(std::ostream &os, const TMatrix &m) {
    for(size_t i = 0; i < m.size_.number_of_rows; ++i) {
        for(size_t j = 0; j < m.size_.number_of_cols; ++j) {
            os << m.matrix_[i][j];
            if(j != m.size_.number_of_cols - 1) {
                os << ' ';
            }
        }
        if(i != m.size_.number_of_rows - 1) {
            os << '\n';
        }
    }
    return os;
}

TMatrix TMatrix::operator+(const TMatrix &other) const {
    assert(this->size_ == other.size_);
    TMatrix result(this->size_);
    for(size_t i = 0; i < size_.number_of_rows; ++i) {
        for(size_t j = 0; j < size_.number_of_cols; ++j) {
            result.matrix_[i][j] = this->matrix_[i][j] + other.matrix_[i][j];
        }
    }
    return result;
}

TMatrix TMatrix::operator*(const TMatrix &other) const {
    assert(AreCompatible(this->size_, other.size_));
    TMatrix result(size_.number_of_rows, other.size_.number_of_cols);
    for(size_t i = 0; i < size_.number_of_rows; ++i) {
        for(size_t j = 0; j < other.size_.number_of_cols; ++j) {
            for(size_t k = 0; k < other.size_.number_of_rows; ++k) {
                result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
            }
        }
    }
    return result;
}

void TMatrix::InitializeMatrix() {
    matrix_.resize(size_.number_of_rows);
    for(size_t i = 0; i < size_.number_of_rows; ++i) {
        matrix_[i].assign(size_.number_of_cols, 0);
    }
}
