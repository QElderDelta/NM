#pragma once

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

#include "matrixops.h"

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
    TMatrix operator-(const TMatrix& other) const;
    TMatrix operator*(const TMatrix& other) const;
    TMatrix operator*(double value) const;

    void Clear();

    double Determinant() const;

    std::vector<double> GetColumn(size_t j) const;

    const std::vector<double>& GetRow(size_t i) const;

    double& GetElement(size_t i, size_t j);

    const double& GetElement(size_t i, size_t j) const;

    double Getl1Norm() const;

    double Getl2Norm() const;

    double GetlinfNorm() const;

    const TMatrixSize& GetSize() const;

    double GetSquaredColumnSum(size_t column, size_t row = 0) const;

    TMatrix InverseMatrix() const;

    std::tuple<TMatrix, TMatrix, TMatrix> LUDecomposition() const;

    std::pair<TMatrix, TMatrix> QRDecomposition() const;

    void SetElement(size_t i, size_t j, double value);

    void SetIdentity();

    void Transpose();
private:
    void SwapCols(size_t i, size_t j);

    void SwapRows(size_t i, size_t j);

    TMatrixSize size_{};
    std::vector<std::vector<double>> matrix_;
};