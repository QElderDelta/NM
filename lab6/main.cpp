#include <hyperbolic_dif_eq_solver.h>

#include <fstream>
#include <iostream>

double left_constraint(double t) {
    return exp(-t) * cos(2 * t);
}

double right_constraint(double t) {
    return 0;
}

double bottom_constraint(double x) {
    return exp(-x) * cos(x);
}

double bottom_constraint_derivative(double x) {
    return -exp(-x) * cos(x);
}

double bottom_constraint_second_derivative(double x) {
    return 2 * sin(x) / exp(x);
}

double free_func(double x, double t) {
    return 0;
}

double analytic_solution(double x, double t) {
    return exp(-t - x) * cos(x) * cos(2 * t);
}

int main() {
    double h, tau, t_max;
    std::cin >> h >> tau >> t_max;
    TaskData t;
    t.left_constraint = left_constraint;
    t.alpha = 0;
    t.beta = 1;
    t.right_constraint = right_constraint;
    t.gamma = 0;
    t.delta = 1;
    t.bottom_constraint = bottom_constraint;
    t.bottom_constraint_derivative = bottom_constraint_derivative;
    t.bottom_constraint_second_derivative = bottom_constraint_second_derivative;
    t.tau = tau;
    t.h = h;
    t.a = 1;
    t.b = 2;
    t.c = -3;
    t.e = 2;
    t.x_left = 0;
    t.x_right = M_PI / 2;
    t.t_bottom = 0;
    t.t_top = t_max;
    t.f = free_func;
    std::ofstream ex("expl.txt");
    auto explicit_solution = solveHyperbolicEquation(t, 0);
    for (int i = 0; i < explicit_solution.GetSize().number_of_rows; ++i) {
        for (int j = 0; j < explicit_solution.GetSize().number_of_cols; ++j) {
            ex << explicit_solution.GetElement(i, j) << ' ';
        }
        ex << '\n';
    }
    std::ofstream im("impl.txt");
    auto implicit_solution = solveHyperbolicEquation(t, 1);
    for (int i = 0; i < implicit_solution.GetSize().number_of_rows; ++i) {
        for (int j = 0; j < implicit_solution.GetSize().number_of_cols; ++j) {
            im << implicit_solution.GetElement(i, j) << ' ';
        }
        im << '\n';
    }
    return 0;
}
