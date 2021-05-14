#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

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

int main() {
    task3_1();
    task3_2();
    return 0;
}
