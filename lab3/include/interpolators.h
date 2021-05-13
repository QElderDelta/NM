#pragma once

#include <functional>
#include <vector>

double lagrangeInterpolation(const std::vector<double>& i_xValues,
                             const std::vector<double>& i_yValues,
                             double i_point);

double newtonInterpolation(const std::vector<double>& i_xValues,
                           const std::vector<double>& i_yValues,
                           double i_point);