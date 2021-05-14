#pragma once

#include <algorithm>
#include <map>
#include <vector>

class DerivativeCalculator {
public:
    //expects sorted vectors
    DerivativeCalculator(const std::vector<double>& i_xValues,
                         const std::vector<double>& i_yValues);
    double getFirstDerivative(double i_point) const;
    double getSecondDerivative(double i_point) const;
private:
    double helper(int index) const;
    std::vector<double> d_xValues;
    std::vector<double> d_yValues;
};