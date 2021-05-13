#include "../include/interpolators.h"

double lagrangeInterpolation(const std::vector<double>& i_xValues,
                             const std::vector<double>& i_yValues,
                             double i_point) {
    double result = 0;
    const double w = [&]() {
        double w = 1;
        for(double i_xValue : i_xValues) {
            w *= (i_point - i_xValue);
        }
        return w;
    }();
    auto getW = [&](int j) {
        double result = 1;
        for(int i = 0; i < i_xValues.size(); ++i) {
            if(i == j) {
                continue;
            }
            result *= (i_xValues[j] - i_xValues[i]);
        }
        return result;
    };
    for(int i = 0; i < i_xValues.size(); ++i) {
        result += i_yValues[i] * w / ((i_point - i_xValues[i]) * getW(i));
    }
    return result;
}

double newtonInterpolation(const std::vector<double> &i_xValues, const std::vector<double> &i_yValues, double i_point) {
    double result = 0;
    std::vector<std::vector<double>> differenceValues(i_xValues.size());
    for(double i_yValue : i_yValues) {
        differenceValues[0].push_back(i_yValue);
    }
    for(int i = 0; i < differenceValues.size() - 1; ++i) {
        for(int j = 0; j < differenceValues[i].size() - 1; ++j) {
            differenceValues[i + 1].push_back((differenceValues[i][j] - differenceValues[i][j + 1])
                                         / (i_xValues[j] - i_xValues[j + i + 1]));
        }
    }
    auto getProduct = [&](int j) {
        double result = 1;
        for(int i = 0; i < j; ++i) {
            result *= (i_point - i_xValues[i]);
        }
        return result;
    };
    for(int i = 0; i < differenceValues.size(); ++i) {
        result += getProduct(i) * differenceValues[i][0];
    }
    return result;
}
