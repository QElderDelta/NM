#include "hyperbolic_dif_eq_solver.h"

namespace {
    double getIthPoint(double start, double h, int i) {
        return start + h * i;
    }

    std::vector<double> singleStep(const TMatrix& res, const TaskData& data, int k, double s) {
        int n = res.GetSize().number_of_cols;
        TTridiagonalMatrix a;
        TMatrix b(n, 1);
        std::stringstream os;
        double a_sqr = pow(data.a, 2);
        double h_sqr = pow(data.h, 2);
        double tau_sqr = pow(data.tau, 2);
        double t_k = getIthPoint(data.t_bottom, data.tau, k - 1);
        double t_k_1 = getIthPoint(data.t_bottom, data.tau, k);
        double value = 0;
        double x_j;
        os << n << '\n';
        os << 2 * a_sqr * data.alpha *
              (1 + h_sqr * (1 + 3 * data.e * data.tau / 2 - data.c * tau_sqr) / (2 * a_sqr * tau_sqr)) / data.h -
              data.beta * (2 * a_sqr - data.h * data.b);
        os << ' ';
        os << -2 * a_sqr * data.alpha / data.h;
        value += 2 * (1 + data.e * data.tau) * res.GetElement(k - 1, 0);
        value -= (1 + data.e * data.tau / 2.) * res.GetElement(k - 2, 0);
        value += tau_sqr * data.f(data.x_left, t_k_1);
        value *= data.alpha * data.h / tau_sqr;
        value -= (2 * a_sqr - data.h * data.b) * data.left_constraint(t_k_1);
        b.SetElement(0, 0, value);
        os << '\n';
        for (int j = 1; j < n - 1; ++j) {
            x_j = getIthPoint(data.x_left, data.h, j);
            os << s * tau_sqr * (-2 * a_sqr + data.b * data.h) / (2 * h_sqr) << ' ';
            os << 1 + 3 * data.e * data.tau / 2. + 2 * a_sqr * tau_sqr * s / h_sqr - data.c * s * tau_sqr << ' ';
            os << s * tau_sqr * (-2 * a_sqr - data.b * data.h) / (2 * h_sqr) << ' ';
            os << '\n';
            value = 0;
            value += (2 + 2 * data.e * data.tau - 2 * a_sqr * tau_sqr * (1 - s) / h_sqr + data.c * (1 - s) * tau_sqr) *
                     res.GetElement(k - 1, j);
            value -= (1 + data.e * data.tau / 2.) * res.GetElement(k - 2, j);
            value += ((1 - s) * tau_sqr * (2 * a_sqr - data.b * data.h) / (2 * h_sqr)) * res.GetElement(k - 1, j - 1);
            value += ((1 - s) * tau_sqr * (2 * a_sqr + data.b * data.h) / (2 * h_sqr)) * res.GetElement(k - 1, j + 1);
            value += (s * data.f(x_j, t_k_1) + (1 - s) * data.f(x_j, t_k)) * tau_sqr;
            b.SetElement(j, 0, value);
        }
        os << -2 * a_sqr * data.gamma / data.h;
        os << ' ';
        os << 2 * a_sqr * data.gamma *
              (1 + h_sqr * (1 + 3 * data.e * data.tau / 2 - data.c * tau_sqr) / (2 * a_sqr * tau_sqr)) / data.h +
              data.delta * (2 * a_sqr + data.h * data.b);
        os << '\n';
        value = 0;
        value += 2 * (1 + data.e * data.tau) * res.GetElement(k - 1, n - 1);
        value -= (1 + data.e * data.tau / 2) * res.GetElement(k - 2, n - 1);
        value += data.f(data.x_right, t_k_1) * tau_sqr;
        value *= data.gamma * data.h / tau_sqr;
        value += (2 * a_sqr + data.h * data.b) * data.right_constraint(t_k_1);
        b.SetElement(n - 1, 0, value);
        os >> a;
        return SweepMethod(a, b).GetColumn(0);
    }
}

TMatrix solveHyperbolicEquation(const TaskData& data, double s) {
    int number_of_rows = (data.t_top - data.t_bottom) / data.tau;
    int number_of_columns = (data.x_right - data.x_left) / data.h + 1;
    TMatrix result(number_of_rows, number_of_columns);
    double x_i;
    for (int i = 0; i < number_of_columns; ++i) {
        x_i = getIthPoint(data.x_left, data.h, i);
        result.SetElement(0, i, data.bottom_constraint(x_i));
        result.SetElement(1, i, data.bottom_constraint(x_i) + data.bottom_constraint_derivative(x_i) * data.tau +
                                pow(data.a, 2) * data.bottom_constraint_second_derivative(x_i) * pow(data.tau, 2) / 2.);
    }
    for (int i = 2; i < number_of_rows; ++i) {
        auto&& row = singleStep(result, data, i, s);
        for (int k = 0; k < number_of_columns; ++k) {
            result.SetElement(i, k, row[k]);
        }
    }
    return result;
}
