#pragma once

#include "../../lab1/lib/matrix.h"
#include "../../lab1/lib/matrixops.h"

#include <cmath>
#include <iostream>
#include <vector>

class LeastSquaresApproximator {
public:
    LeastSquaresApproximator(const std::vector<double>& i_xValues,
                             const std::vector<double>& i_yValues,
                             int order,
                             std::ostream& o_logStream = std::cout);
    double approximate(double i_point) const;
    std::vector<double> getCoefficients() const;
private:
    std::vector<double> d_coefficients;
};