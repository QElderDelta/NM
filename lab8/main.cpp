#include <2d_eq_solver.h>

#include <fstream>
#include <iostream>

double leftConstraint(double y, double t) {
    return 0;
}

double rightConstraint(double y, double t) {
    return y * cos(t);
}

double bottomConstraint(double x, double t) {
    return 0;
}

double topConstraint(double x, double t) {
    return x * cos(t);
}

double leftEqPart(double x, double y, double t) {
    return -x * y * sin(t);
}

double firstLayer(double x, double y) {
    return x * y;
}

double getIthPoint(double start, double h, int i) {
    return start + h * i;
}

double analytic(double x, double y, double t) {
    return x * y * cos(t);
}

int main() {
    double hx, hy, tau, t_end;

    std::cin >> hx >> hy >> tau >> t_end;

    TaskData data;
    data.left_constraint = leftConstraint;
    data.right_constraint = rightConstraint;
    data.bottom_constraint = bottomConstraint;
    data.top_constraint = topConstraint;
    data.left_part = leftEqPart;
    data.first_layer_constraint = firstLayer;
    data.x_left = 0;
    data.x_right = 1;
    data.y_bottom = 0;
    data.y_top = 1;
    data.t_bottom = 0;
    data.t_top = t_end;
    data.hx = hx;
    data.hy = hy;
    data.tau = tau;

    auto var_res = varyingDirectionsMethod(data);
//    int n = (data.y_top - data.y_bottom) / data.hy;
//    int m = (data.x_right - data.x_left) / data.hx;
//    std::cout << res[0] << '\n';
//    std::cout << "----------------------" << '\n';
//    for (int i = 0; i < n + 1; ++i) {
//        for (int j = 0; j < m + 1; ++j) {
//            std::cout << analytic(getIthPoint(data.x_left, data.hx, j), getIthPoint(data.y_bottom, data.hy, i), 0)
//                      << ' ';
//        }
//        std::cout << '\n';
//    }
//    std::cout << "-----------------" << '\n';
//    std::cout << res[1] << '\n';
//    std::cout << "----------------------" << '\n';
//    for (int i = 0; i < n + 1; ++i) {
//        for (int j = 0; j < m + 1; ++j) {
//            auto xd = analytic(getIthPoint(data.x_left, data.hx, j), getIthPoint(data.y_bottom, data.hy, i), 0.5);
//            std::cout << analytic(getIthPoint(data.x_left, data.hx, j), getIthPoint(data.y_bottom, data.hy, i), 0.5)
//                      << ' ';
//        }
//        std::cout << '\n';
//    }
    std::ofstream var("var.txt");
    for (const auto& item: var_res) {
        var << item << '\n';
    }
    auto frac_res = fractionalStepsMethod(data);
    std::ofstream frac("frac.txt");
    for (const auto& item: frac_res) {
        frac << item << '\n';
    }
    return 0;
}
