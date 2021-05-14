#include <cassert>

#include "../include/derivative.h"


DerivativeCalculator::DerivativeCalculator(const std::vector<double> &i_xValues,
                                           const std::vector<double> &i_yValues) : d_xValues(i_xValues),
                                                                                   d_yValues(i_yValues) {
    assert(i_xValues.size() == i_yValues.size());
}

double DerivativeCalculator::getFirstDerivative(double i_point) const {
    if(i_point <= d_xValues.front() || i_point >= d_xValues.back()) {
        throw std::invalid_argument("Point is out of range");
    }
    auto it = std::lower_bound(d_xValues.begin(), d_xValues.end(), i_point);
    int index = std::distance(d_xValues.begin(), it);
    if(i_point == *it) {
        auto first = helper(index - 1);
        return first + (helper(index) - first) * (d_xValues[index] - d_xValues[index - 1])
                     / (d_xValues[index + 1] + d_xValues[index - 1]);
    } else {
        return helper(index);
    }
}

double DerivativeCalculator::getSecondDerivative(double i_point) const {
    if(i_point <= d_xValues.front() || i_point >= d_xValues.back()) {
        throw std::invalid_argument("Point is out of range");
    }
    auto it = std::lower_bound(d_xValues.begin(), d_xValues.end(), i_point);
    int index = std::distance(d_xValues.begin(), it);
    if(i_point == *it) {
        return 2 * (helper(index) - helper(index - 1)) / (d_xValues[index + 1] - d_xValues[index - 1]);
    } else {
        if(d_xValues.size() - 1 - index < 1) {
            throw std::invalid_argument("Can't calculate second derivative for this point");
        }
        return 2 * (helper((index + 1) + helper(index))) / (d_xValues[index + 2] - d_xValues[index]);
    }
}

double DerivativeCalculator::helper(int index) const {
    return (d_yValues[index + 1] - d_yValues[index]) / (d_xValues[index + 1] - d_xValues[index]);
}
