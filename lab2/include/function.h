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


class ITwoArgumentFunction {
public:
    virtual double getValue(double i_x1, double i_x2) = 0;
    virtual double getFirstDerivativeByFirstArgument(double i_x1, double i_x2) = 0;
    virtual double getFirstDerivativeBySecondArgument(double i_x1, double i_x2) = 0;
    ~ITwoArgumentFunction() = default;
};

class SecondTaskFunction1 : public ITwoArgumentFunction {
public:
    double getValue(double i_x1, double i_x2) override;
    double getFirstDerivativeByFirstArgument(double i_x1, double i_x2) override;
    double getFirstDerivativeBySecondArgument(double i_x1, double i_x2) override;
};

class SecondTaskFunction2 : public ITwoArgumentFunction {
public:
    double getValue(double i_x1, double i_x2) override;
    double getFirstDerivativeByFirstArgument(double i_x1, double i_x2) override;
    double getFirstDerivativeBySecondArgument(double i_x1, double i_x2) override;
};

class SecondTaskFunction1Transformed : public SecondTaskFunction1 {
    double getValue(double i_x1, double i_x2) override;
};

class SecondTaskFunction2Transformed : public SecondTaskFunction2 {
    double getValue(double i_x1, double i_x2) override;
};