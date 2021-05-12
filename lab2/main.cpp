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

int main() {
    task2_1();
    return 0;
}
