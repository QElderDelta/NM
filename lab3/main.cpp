#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "include/approximators.h"
#include "include/derivative.h"
#include "include/integral.h"
#include "include/interpolators.h"

void task3_1() {
    std::ifstream is("../task3_1.txt");
    std::ofstream os("../task3_1_result.txt");
    int n;
    double point;
    is >> n >> point;
    std::vector<double> points1(n);
    std::vector<double> values1(n);
    std::vector<double> points2(n);
    std::vector<double> values2(n);
    for(int i = 0; i < n; ++i) {
        is >> points1[i];
        values1[i] = log(points1[i]);
    }
    for(int i = 0; i < n; ++i) {
        is >> points2[i];
        values2[i] = log(points2[i]);
    }
    os << "Lagrange's interpolation result for the first set of points: ";
    os << lagrangeInterpolation(points1, values1, point) << '\n';
    os << "Actual value is " << log(point) << '\n';
    os << "Absolute error is " << std::abs(log(point) - lagrangeInterpolation(points1, values1, point)) << '\n';
    os << "Newton's interpolation result for the second set of points: ";
    os << newtonInterpolation(points2, values2, point) << '\n';
    os << "Actual value is " << log(point) << '\n';
    os << "Absolute error is " << std::abs(log(point) - newtonInterpolation(points2, values2, point)) << '\n';
}

void task3_2() {
    std::ifstream is("../task3_2.txt");
    std::ofstream os("../task3_2_result.txt");
    std::ofstream plotData("../task3_2_plot.txt");
    int n;
    double point;
    is >> n >> point;
    std::vector<double> points(n);
    std::vector<double> values(n);
    for(int i = 0; i < n; ++i) {
        is >> points[i];
    }
    for(int i = 0; i < n; ++i) {
        is >> values[i];
    }
    SplineInterpolator s(points, values);
    os << "Result of interpolation: ";
    os << s.getValue(point) << '\n';
    s.getSplineInfo(os);
    s.getSplineInfoWithoutText(plotData);
}

void task3_3() {
    std::ifstream is("../task3_3.txt");
    std::ofstream os("../task3_3_result.txt");
    std::ofstream plotData("../task3_3_plot.txt");
    int n;
    is >> n;
    std::vector<double> points(n);
    std::vector<double> values(n);
    for(int i = 0; i < n; ++i) {
        is >> points[i];
    }
    for(int i = 0; i < n; ++i) {
        is >> values[i];
    }
    os << "Approximation using first order polynom" << '\n';
    LeastSquaresApproximator firstOrder(points, values, 1, os);
    os << "Resulting function: ";
    auto aCoeffs = firstOrder.getCoefficients();
    plotData << aCoeffs[0] << " " << aCoeffs[1] << '\n';
    os << aCoeffs[0] << " + " << aCoeffs[1] << " * x" << '\n';
    os << "Approximation using second order polynom" << '\n';
    LeastSquaresApproximator secondOrder(points, values, 2, os);
    os << "Resulting function: ";
    auto bCoeffs = secondOrder.getCoefficients();
    plotData << bCoeffs[0] << " " << bCoeffs[1] << " " << bCoeffs[2] << '\n';
    os << bCoeffs[0] << " + " << bCoeffs[1] << " * x + " << bCoeffs[2] << " * x^2" << '\n';
}

void task3_4() {
    std::ifstream is("../task3_4.txt");
    std::ofstream os("../task3_4_result.txt");
    int n;
    double point;
    is >> n >> point;
    std::vector<double> points(n);
    std::vector<double> values(n);
    for(int i = 0; i < n; ++i) {
        is >> points[i];
    }
    for(int i = 0; i < n; ++i) {
        is >> values[i];
    }
    DerivativeCalculator c(points, values);
    os << "First derivative: " << c.getFirstDerivative(point) << '\n';
    os << "Second derivative: " << c.getSecondDerivative(point) << '\n';
}

void task3_5() {
    std::ifstream is("../task3_5.txt");
    std::ofstream os("../task3_5_result.txt");
    double start, end, h1, h2;
    is >> start >> end >> h1 >> h2;
    auto f = [](double x) {
        return 1 / ((2 * x + 7) * (3 * x + 4));
    };
    const double exactValue = 0.104471;
    os << "Exact value is " << exactValue << '\n';
    os << "Result of integration using rectangle method with step " << h1 << ": ";
    os << integrateUsingRectangleMethod(f, start, end, h1) << '\n';
    os << "Result of integration using rectangle method with step " << h2 << ": ";
    os << integrateUsingRectangleMethod(f, start, end, h2) << '\n';
    os << "Result of integration using trapezoid method with step " << h1 << ": ";
    os << integrateUsingTrapezoidMethod(f, start, end, h1) << '\n';
    os << "Result of integration using trapezoid method with step " << h2 << ": ";
    os << integrateUsingTrapezoidMethod(f, start, end, h2) << '\n';
    os << "Result of integration using Simpson's method with step " << h1 << ": ";
    os << integrateUsingSimpsonMethod(f, start, end, h1) << '\n';
    os << "Result of integration using Simpson's method with step " << h2 << ": ";
    os << integrateUsingSimpsonMethod(f, start, end, h2) << '\n';
    os << "Error for rectangle method counted using Runge-Romberg method: ";
    auto res1 = rungeRombergMethod(integrateUsingRectangleMethod(f, start, end, h1),
                                  integrateUsingRectangleMethod(f, start, end, h2),
                                  h1 / h2,
                                  2);
    os << res1 << '\n';
    auto res2 = rungeRombergMethod(integrateUsingTrapezoidMethod(f, start, end, h1),
                                   integrateUsingTrapezoidMethod(f, start, end, h2),
                                   h1 / h2,
                                   2);
    os << "Error for trapezoid method counted using Runge-Romberg method: ";
    os << res2 << '\n';
    auto res3 = rungeRombergMethod(integrateUsingSimpsonMethod(f, start, end, h1),
                                   integrateUsingSimpsonMethod(f, start, end, h2),
                                   h1 / h2,
                                   4);
    os << "Error for Simpson's method counted using Runge-Romberg method: ";
    os << res3 << '\n';
}

int main() {
    task3_1();
    task3_2();
    task3_3();
    task3_4();
    task3_5();
    return 0;
}
