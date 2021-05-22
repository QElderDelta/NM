#pragma once

#include "../../lab2/include/function.h"

#include <vector>

class Function {
public:
    Function(int order);
    int getOrder() const;
    virtual double getValue(const std::vector<double>& i_values) const = 0;
private:
    int d_order;
};

struct Point {
    double x, y, z;

    Point(double x, double y, double z) : x(x), y(y), z(z) {};
};

std::vector<std::pair<double, double>> eulerMethod(ITwoArgumentFunction& i_f,
                                                   double i_a,
                                                   double i_b,
                                                   double i_h,
                                                   double i_leftConstraint);

std::vector<std::pair<double, double>> rungeKuttaFourthOrderMethod(ITwoArgumentFunction& i_f,
                                                                   double i_a,
                                                                   double i_b,
                                                                   double i_h,
                                                                   double i_leftConstraint);

std::vector<std::pair<double, double>> adamsFourthOrderMethod(ITwoArgumentFunction& i_f,
                                                              double i_a,
                                                              double i_b,
                                                              double i_h,
                                                              double i_leftConstraint);

std::vector<Point> solveSecondOrderUsingEulerMethod(Function& i_f,
                                                    double i_a,
                                                    double i_b,
                                                    double i_h,
                                                    double i_yLeftConstraint,
                                                    double i_zLeftConstraint);

std::vector<Point> solveSecondOrderUsingRungeKuttaFourthOrderMethod(Function& i_f,
                                                                    double i_a,
                                                                    double i_b,
                                                                    double i_h,
                                                                    double i_yLeftConstraint,
                                                                    double i_zLeftConstraint);

std::vector<Point> solveSecondOrderUsingAdamsFourthOrderMethod(Function& i_f,
                                                               double i_a,
                                                               double i_b,
                                                               double i_h,
                                                               double i_yLeftConstraint,
                                                               double i_zLeftConstraint);
