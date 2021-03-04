#pragma once

#include <cassert>
#include <iostream>
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
private:
    void InitializeMatrix();

    TMatrixSize size_;
    std::vector<std::vector<double>> matrix_;
};