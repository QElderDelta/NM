#include <fstream>
#include <iostream>

#include "include/function.h"
#include "include/methods.h"

void task2_1() {
    std::ifstream is("../task2_1.txt");
    std::ofstream os("../task2_1_result.txt");
    double eps, x0, q;
    is >> eps;
    is >> x0;
    TaskFunction f;
    TransformedFunction tf;
    os << "Root found by Newton's method:" << '\n';
    os << newtonsMethod(f, x0, eps, os) << '\n';
    is >> x0;
    is >> q;
    os << "Root found by simple iterations method:" << '\n';
    os << simpleIterationsMethod(tf, x0, q, eps, os) << '\n';
}

void task2_2() {
    std::ifstream is("../task2_2.txt");
    std::ofstream os("../task2_2_result.txt");
    double eps, x0_1, x0_2, q;
    is >> eps;
    is >> x0_1 >> x0_2;
    SecondTaskFunction1 f1;
    SecondTaskFunction2 f2;
    os << "Solution found by Newton's method: " << '\n';
    auto [x1, x2] = newtonsMethod(f1, f2, {x0_1, x0_2}, eps, os);
    os << "x1: " << x1 << " x2: " << x2 << '\n';
    is >> q;
    os << "Solution found by simple iterations method: " << '\n';
    SecondTaskFunction1Transformed f3;
    SecondTaskFunction2Transformed f4;
    auto [x3, x4] = simpleIterationsMethod(f3, f4, {x0_1, x0_2}, q, eps, os);
    os << "x1: " << x3 << " x2: " << x4 << '\n';
}

int main() {
    task2_1();
    task2_2();
    return 0;
}
