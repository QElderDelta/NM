#pragma once

#include "../../lab1/lib/matrixops.h"

#include <map>
#include <sstream>
#include <vector>

double lagrangeInterpolation(const std::vector<double>& i_xValues,
                             const std::vector<double>& i_yValues,
                             double i_point);

double newtonInterpolation(const std::vector<double>& i_xValues,
                           const std::vector<double>& i_yValues,
                           double i_point);

class SplineInterpolator {
public:
    SplineInterpolator(const std::vector<double>& i_xValues, const std::vector<double>& i_yValues);
    double getValue(double i_point) const;
    void getSplineInfo(std::ostream& os) const;
private:
    std::map<double, std::vector<double>> d_splines;
};