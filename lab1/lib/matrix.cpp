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
    SetIdentity();
}

TMatrix::TMatrix(const TMatrixSize &size) : size_(size) {
    SetIdentity();
}

std::istream &operator>>(std::istream &is, TMatrix &m) {
    is >> m.size_.number_of_rows >> m.size_.number_of_cols;
    m.SetIdentity();
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
    if(this->size_ != other.size_) {
        throw std::invalid_argument("Matrices don't have same sizes");
    }
    TMatrix result(this->size_);
    for(size_t i = 0; i < size_.number_of_rows; ++i) {
        for(size_t j = 0; j < size_.number_of_cols; ++j) {
            result.matrix_[i][j] = this->matrix_[i][j] + other.matrix_[i][j];
        }
    }
    return result;
}

TMatrix TMatrix::operator*(const TMatrix &other) const {
    if(!AreCompatible(this->size_, other.size_)) {
        throw std::invalid_argument("Matrices aren't compatible");
    }
    TMatrix result(size_.number_of_rows, other.size_.number_of_cols);
    result.Clear();
    for(size_t i = 0; i < size_.number_of_rows; ++i) {
        for(size_t j = 0; j < other.size_.number_of_cols; ++j) {
            for(size_t k = 0; k < other.size_.number_of_rows; ++k) {
                result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
            }
        }
    }
    return result;
}

void TMatrix::Clear() {
    matrix_.resize(size_.number_of_rows);
    for(size_t i = 0; i < size_.number_of_rows; ++i) {
        matrix_[i].assign(size_.number_of_cols, 0);
    }
}

void TMatrix::SetIdentity() {
    matrix_.resize(size_.number_of_rows);
    for(size_t i = 0; i < size_.number_of_rows; ++i) {
        matrix_[i].resize(size_.number_of_cols);
        for(int j = 0; j < size_.number_of_cols; ++j) {
            if(i == j) {
                matrix_[i][j] = 1;
            } else {
                matrix_[i][j] = 0;
            }
        }
    }
}

std::tuple<TMatrix, TMatrix, TMatrix> TMatrix::LUDecomposition() const {
    TMatrix l(size_);
    TMatrix u(*this);
    TMatrix m(size_);
    TMatrix p(size_);
    std::vector<std::pair<size_t, size_t>> swaps;
    std::vector<size_t> rows_pos(size_.number_of_rows);
    std::iota(rows_pos.begin(), rows_pos.end(), 0);
    for(size_t k = 0; k < size_.number_of_cols - 1; ++k) {
        double max_element = std::abs(u.matrix_[k][k]);
        size_t max_element_row = k;
        for(size_t i = k + 1; i < size_.number_of_rows; ++i) {
            if(std::abs(u.matrix_[i][k]) > max_element) {
                max_element = std::abs(u.matrix_[i][k]);
                max_element_row = i;
            }
        }
        if(max_element_row != k) {
            std::swap(u.matrix_[k], u.matrix_[max_element_row]);
            std::swap(rows_pos[k], rows_pos[max_element_row]);
            l.SwapRows(k, max_element_row);
            l.SwapCols(k, max_element_row);
            swaps.emplace_back(k, max_element_row);
        }
        for(size_t i = k + 1; i < size_.number_of_rows; ++i) {
            m.matrix_[i][k] = -u.matrix_[i][k] / u.matrix_[k][k];
        }
        u = m * u;
        for(size_t i = k + 1; i < size_.number_of_rows; ++i) {
            l.matrix_[i][k] = m.matrix_[i][k] * (-1);
        }
        m.SetIdentity();
    }
    for(auto it = swaps.rbegin(); it != swaps.rend(); ++it) {
        p.SwapRows(it->first, it->second);
    }
    return {l, u, p};
}

TMatrix TMatrix::InverseMatrix() const {
    return TMatrix();
}

double TMatrix::Determinant() const {
    return 0;
}

void TMatrix::SwapRows(size_t i, size_t j) {
    if(i >= size_.number_of_rows || j >= size_.number_of_cols) {
        throw std::out_of_range("");
    }
    std::swap(matrix_[i], matrix_[j]);
}

void TMatrix::SwapCols(size_t i, size_t j) {
    if(i >= size_.number_of_rows || j >= size_.number_of_cols) {
        throw std::out_of_range("");
    }
    for(size_t k = 0; k < size_.number_of_rows; ++k) {
        std::swap(matrix_[k][i], matrix_[k][j]);
    }
}

double &TMatrix::GetElement(size_t i, size_t j) {
    if(i >= size_.number_of_rows || j >= size_.number_of_cols) {
        throw std::out_of_range("");
    }
    return matrix_[i][j];
}

const double &TMatrix::GetElement(size_t i, size_t j) const {
    if(i >= size_.number_of_rows || j >= size_.number_of_cols) {
        throw std::out_of_range("");
    }
    return matrix_[i][j];
}

void TMatrix::SetElement(size_t i, size_t j, double value) {
    if(i >= size_.number_of_rows || j >= size_.number_of_cols) {
        throw std::out_of_range("");
    }
    matrix_[i][j] = value;
}
