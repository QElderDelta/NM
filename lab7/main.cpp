#include <elliptic_eq_solver.h>

#include <fstream>
#include <iostream>

double leftConstraint(double y) {
    return cos(y);
}

double rightConstraint(double y) {
    return 0;
}

double bottomConstraint(double x) {
    return cos(x);
}

double topConstraint(double x) {
    return 0;
}

double leftEqPart(double u) {
    return -2 * u;
}

double getIthPoint(double start, double h, int i) {
    return start + h * i;
}

int main() {
    double hx, hy, eps;

    std::cin >> hx >> hy >> eps;

    TaskData data;
    data.left_constraint = leftConstraint;
    data.right_constraint = rightConstraint;
    data.bottom_constraint = bottomConstraint;
    data.top_constraint = topConstraint;
    data.left_eq_part = leftEqPart;
    data.x_left = 0;
    data.x_right = M_PI / 2;
    data.y_bottom = 0;
    data.y_top = M_PI / 2;
    data.eps = eps;
    data.hx = hx;
    data.hy = hy;

    std::ofstream lieb("lieb.txt");
    std::ofstream lieb_eps("lieb_eps.txt");
    auto lieb_solution = liebmannMethod(data, &lieb_eps);
    lieb << lieb_solution << '\n';
    std::ofstream seid("seid.txt");
    std::ofstream seid_eps("seid_eps.txt");
    auto seid_solution = seidelMethod(data, &seid_eps);
    seid << seid_solution << '\n';
    std::ofstream relax("relax.txt");
    std::ofstream relax_eps("relax_eps.txt");
    auto relax_solution = relaxationMethod(data, 0.9, &relax_eps);
    relax << relax_solution << '\n';
    return 0;
}
