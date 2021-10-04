#pragma once

#include <iostream>
#include <vector>

class TTridiagonalMatrix {
public:
    TTridiagonalMatrix() = default;
    friend std::istream& operator>>(std::istream& is, TTridiagonalMatrix& m);
    double& GetElement(size_t i, size_t j);
    const double& GetElement(size_t i, size_t j) const;
    void SetElement(size_t i, size_t j, double value);
private:
    size_t number_of_rows_{};
    std::vector<std::vector<double>> matrix_;
};
