#include <matrix.h>
#include <parabolic_dif_eq_solver.h>

#include <fstream>
#include <iostream>

double left_constraint(double t) {
    return sin(t);
}

double right_constraint(double t) {
    return -sin(t);
}

double bottom_constraint(double x) {
    return 0;
}

double free_func(double x, double t) {
    return cos(x) * (cos(t) + sin(t));
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        throw std::runtime_error("Grid parameters weren't passed to the function");
    }
    double h = std::stod(argv[1]);
    double tau = std::stod(argv[2]);
    std::ofstream f1("expl.txt", std::ios::trunc);
    std::ofstream f2("impl.txt", std::ios::trunc);
    std::ofstream f3("exim.txt", std::ios::trunc);
    std::ofstream f4("analytic.txt", std::ios::trunc);
    TaskData d{
            left_constraint,
            0, //alpha
            1, //beta
            right_constraint,
            1, //gamma
            0, //delta
            bottom_constraint,
            tau, //tau,
            h, //h
            1, //a
            0, //x_left
            M_PI / 2, //x_right
            0, //t_bottom
            10, //t_top
            free_func
    };
    auto ex_res = CrankNicolsonSolution(d, 0);
    f1 << ex_res << '\n';
    auto imp_res = CrankNicolsonSolution(d, 1);
    f2 << imp_res << '\n';
    auto crank_nicl_res = CrankNicolsonSolution(d, 0.5);
    f3 << crank_nicl_res << '\n';
    for(int i = 0; i < (d.t_top - d.t_bottom) / d.tau; ++i) {
        for(int j = 0; j < (d.x_right - d.x_left) / d.h + 1; ++j) {
            f4 << sin(GetIthPoint(d.t_bottom, d.tau, i)) * cos(GetIthPoint(d.x_left, d.h, j)) << ' ';
        }
        f4 << '\n';
    }
    return 0;
}
