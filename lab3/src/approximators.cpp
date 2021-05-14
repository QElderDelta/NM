#include "../include/approximators.h"

LeastSquaresApproximator::LeastSquaresApproximator(const std::vector<double> &i_xValues,
                                                   const std::vector<double> &i_yValues, int order,
                                                   std::ostream &o_logStream) {
    assert(order > 0);
    auto getXSum = [&](int power) {
        double result = 0;
        for(const auto& x : i_xValues) {
            result += pow(x, power);
        }
        return result;
    };
    auto getYSum = [&](int power) {
        double result = 0;
        for(int i = 0; i < i_xValues.size(); ++i) {
            result += i_yValues[i] * pow(i_xValues[i], power);
        }
        return result;
    };
    TMatrix a(order + 1, order + 1);
    TMatrix b(order + 1, 1);
    for(int i = 0; i < order + 1; ++i) {
        for(int j = 0; j < order + 1; ++j) {
            a.SetElement(i, j, getXSum(i + j));
        }
    }
    for(int i = 0; i < order + 1; ++i) {
        b.SetElement(i, 0, getYSum(i));
    }
    TMatrix c = SolveLinearSystem(a, b);
    for(int i = 0; i < order + 1; ++i) {
        d_coefficients.push_back(c.GetElement(i, 0));
    }
    double sumOfResidualsSquared = 0;
    for(int i = 0; i < i_xValues.size(); ++i) {
        sumOfResidualsSquared += pow(i_yValues[i] - approximate(i_xValues[i]), 2);
    }
    o_logStream << "Sum of the squares of residuals is " << sumOfResidualsSquared << '\n';
}

double LeastSquaresApproximator::approximate(double i_point) const {
    double result = 0;
    for(int i = 0; i < d_coefficients.size(); ++i) {
        result += d_coefficients[i] * pow(i_point, i);
    }
    return result;
}

std::vector<double> LeastSquaresApproximator::getCoefficients() const {
    return d_coefficients;
}
