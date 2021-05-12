#pragma once

#include <cmath>

class ISingleArgumentFunction {
public:
    virtual double getValue(double i_x) = 0;
    virtual double firstDerivative(double i_x) = 0;
    ~ISingleArgumentFunction() = default;
};

class TaskFunction : public ISingleArgumentFunction {
    double getValue(double i_x) override;
    double firstDerivative(double i_x) override;
};

class TransformedFunction : public TaskFunction {
    double getValue(double i_x) override;
};
