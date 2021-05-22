#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "include/cauchyproblemsolvers.h"
#include "../lab3/include/integral.h"

class FirstTaskFunction : public Function {
    using Function::Function;

    double getValue(const std::vector<double>& i_values) const override {
        assert(i_values.size() == this->getOrder());
        return (1 + 2 * pow(tan(i_values[0]), 2)) * i_values[1];
    }
};



void task4_1() {
    std::ofstream os("../task4_1_result.txt");
    FirstTaskFunction f(3);
    const auto exactSolution = []() {
        std::vector<std::pair<double, double>> res;
        double current = 0;
        auto f = [](double x) {
            return 1 / cos(x) + sin(x) + x / cos(x);
        };
        while(current <= 1) {
            res.emplace_back(current, f(current));
            current += 0.1;
        }
        return res;
    }();
    os << std::fixed << std::setprecision(7);
    os << "Exact solution:" << '\n';
    for(const auto& item : exactSolution) {
        os << "{x: " << item.first << ", y: " << item.second << "}" << '\n';
    }
    os << "Solution using Euler's method:" << '\n';
    auto res = solveSecondOrderUsingEulerMethod(f, 0, 1, 0.1, 1, 2);
    for(const auto& item : res) {
        os << "{x: " << item.x << ", y: " << item.y << "}" << '\n';
    }
    os << "In point x = 1 absolute error for Euler's method is: ";
    os << std::abs(exactSolution.back().second - res.back().y) << '\n';
    auto tmp = res.back().y;
    os << "In point x = 1 error for Euler's method counted using Runge-Romberg method: ";
    res = solveSecondOrderUsingEulerMethod(f, 0, 1, 0.05, 1, 2);
    os << rungeRombergMethod(tmp, res.back().y, 2, 1) << '\n';
    os << "Solution using Runge-Kutta fourth order method:" << '\n';
    res = solveSecondOrderUsingRungeKuttaFourthOrderMethod(f, 0, 1, 0.1, 1, 2);
    for(const auto& item : res) {
        os << "{x: " << item.x << ", y: " << item.y << "}" << '\n';
    }
    os << "In point x = 1 absolute error for Runge-Kutta fourth order method is: ";
    os << std::abs(exactSolution.back().second - res.back().y) << '\n';
    tmp = res.back().y;
    os << "In point x = 1 error for Runge-Kutta fourth order method counted using Runge-Romberg method: ";
    res = solveSecondOrderUsingRungeKuttaFourthOrderMethod(f, 0, 1, 0.05, 1, 2);
    os << rungeRombergMethod(tmp, res.back().y, 2, 4) << '\n';
    os << "Solution using Adams' fourth order method:" << '\n';
    res = solveSecondOrderUsingAdamsFourthOrderMethod(f, 0, 1, 0.1, 1, 2);
    for(const auto& item : res) {
        os << "{x: " << item.x << ", y: " << item.y << "}" << '\n';
    }
    os << "In point x = 1 absolute error for Adams' fourth order method is: ";
    os << std::abs(exactSolution.back().second - res.back().y) << '\n';
    tmp = res.back().y;
    os << "In point x = 1 error for Adams' fourth order method counted using Runge-Romberg method: ";
    res = solveSecondOrderUsingAdamsFourthOrderMethod(f, 0, 1, 0.05, 1, 2);
    os << rungeRombergMethod(tmp, res.back().y, 2, 4) << '\n';
}

int main() {
    task4_1();
    return 0;
}
