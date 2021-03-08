#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

struct TMatrixSize {
    size_t number_of_rows, number_of_cols;

    bool operator==(const TMatrixSize& other) const;
    bool operator!=(const TMatrixSize& other) const;
};

bool AreCompatible(const TMatrixSize& first, const TMatrixSize& second);

class TMatrix {
public:
    TMatrix() = default;
    TMatrix(size_t number_of_rows, size_t number_of_cols);
    explicit TMatrix(const TMatrixSize& size);

    friend std::istream& operator>>(std::istream& is, TMatrix& m);
    friend std::ostream& operator<<(std::ostream& is, const TMatrix& m);

    TMatrix operator+(const TMatrix& other) const;
    TMatrix operator*(const TMatrix& other) const;

    std::tuple<TMatrix, TMatrix, TMatrix> LUDecomposition() const;

    TMatrix InverseMatrix() const;

    double Determinant() const;

    void Clear();

    void SetIdentity();

    double& GetElement(size_t i, size_t j);

    const double& GetElement(size_t i, size_t j) const;

    void SetElement(size_t i, size_t j, double value);
private:
    void SwapRows(size_t i, size_t j);

    void SwapCols(size_t i, size_t j);

    TMatrixSize size_{};
    std::vector<std::vector<double>> matrix_;
};