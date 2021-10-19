#include "parabolic_dif_eq_solver.h"

double GetIthPoint(double start, double h, int i) {
    return start + h * i;
}

std::vector<double>
ExplicitSchemaSolutionSingleStep(const TMatrix& res, const TaskData& data, int k, std::optional<double> theta) {
    int n = res.GetSize().number_of_cols;
    std::vector<double> result(n);
    double sigma = data.a * data.a * data.tau / data.h / data.h;
    double t_k = GetIthPoint(data.t_bottom, data.tau, k);
    for (int j = 1; j < n - 1; ++j) {
        auto current_x = GetIthPoint(data.x_left, data.h, j);
        result[j] =
                sigma * res.GetElement(k - 1, j + 1) + (1 - 2 * sigma) * res.GetElement(k - 1, j) +
                sigma * res.GetElement(k - 1, j - 1) +
                data.tau * data.f(current_x, t_k);
    }
    result[0] = -(data.alpha / data.h) * result[1] / (data.beta - data.alpha / data.h);
    result[0] += data.left_constraint(t_k) / (data.beta - data.alpha / data.h);
    result[n - 1] = (data.right_constraint(t_k) * data.h + data.gamma * result[n - 2]) /
                    (data.gamma + data.h * data.delta);
    if (theta) {
        for (int i = 0; i < n; ++i) {
            result[i] *= theta.value();
        }
    }
    return result;
}

TMatrix ExplicitSchemaSolution(const TaskData& data) {
    double sigma = data.a * data.a * data.tau / data.h / data.h;
    if (sigma > 0.5) {
        throw std::invalid_argument("Sigma is more that 1/2, method is unstable");
    }
    int number_of_rows = (data.t_top - data.t_bottom) / data.tau;
    int number_of_columns = (data.x_right - data.x_left) / data.h + 1;
    TMatrix result(number_of_rows, number_of_columns);
    for (int i = 0; i < number_of_columns; ++i) {
        result.SetElement(0, i, data.bottom_constraint(GetIthPoint(data.x_left, data.h, i)));
    }
    for (int i = 1; i < number_of_rows; ++i) {
        auto row = ExplicitSchemaSolutionSingleStep(result, data, i);
        for (int k = 0; k < number_of_columns; ++k) {
            result.SetElement(i, k, row[k]);
        }
    }
    return result;
}

std::vector<double> ImplicitSchemaSolutionSingleStep(const TMatrix& res, const TaskData& data, int k,
                                                     std::optional<std::vector<double>> to_add,
                                                     std::optional<double> theta) {
    int n = res.GetSize().number_of_cols;
    double sigma = data.a * data.a * data.tau / data.h / data.h;
    std::vector<double> result;
    TTridiagonalMatrix a;
    TMatrix b(n, 1);
    std::stringstream os;
    double t_k = GetIthPoint(data.t_bottom, data.tau, k);
    os << n << '\n';
    os << data.beta - data.alpha * data.h << ' ' << data.alpha / data.h << '\n';
    b.SetElement(0, 0, data.left_constraint(t_k));
    for (int j = 1; j < n - 1; ++j) {
        auto current_x = GetIthPoint(data.x_left, data.h, j);
        os << sigma << ' ' << -(1 + 2 * sigma) << ' ' << sigma << '\n';
        if (to_add) {
            b.SetElement(j, 0, -(res.GetElement(k - 1, j) + data.tau * data.f(current_x, t_k) + to_add.value()[j]));
        } else {
            b.SetElement(j, 0, -(res.GetElement(k - 1, j) + data.tau * data.f(current_x, t_k)));
        }
    }
    os << -data.gamma / data.h << ' ' << data.delta + data.gamma / data.h << '\n';
    b.SetElement(n - 1, 0, data.right_constraint(t_k));
    os >> a;
    auto c = SweepMethod(a, b);
    result = c.GetColumn(0);
    if (theta) {
        for (int i = 0; i < n; ++i) {
            result[i] *= theta.value();
        }
    }
    return result;
}

TMatrix ImplicitSchemaSolution(const TaskData& data) {
    int number_of_rows = (data.t_top - data.t_bottom) / data.tau;
    int number_of_columns = (data.x_right - data.x_left) / data.h + 1;
    TMatrix result(number_of_rows, number_of_columns);
    for (int i = 0; i < number_of_columns; ++i) {
        result.SetElement(0, i, data.bottom_constraint(GetIthPoint(data.x_left, data.h, i)));
    }
    for (int i = 1; i < number_of_rows; ++i) {
        auto row = ImplicitSchemaSolutionSingleStep(result, data, i);
        for (int k = 0; k < number_of_columns; ++k) {
            result.SetElement(i, k, row[k]);
        }
    }
    return result;
}

std::vector<double> CrankNicolsonSolutionSingleStep(const TMatrix& res, const TaskData& data, int k, double theta) {
    int n = res.GetSize().number_of_cols;
    double sigma = data.a * data.a * data.tau / data.h / data.h;
    std::vector<double> result;
    TTridiagonalMatrix a;
    TMatrix b(n, 1);
    std::stringstream os;
    double t_k = GetIthPoint(data.t_bottom, data.tau, k);
    auto getLowerRowSum = [&res, k](int j) {
        return res.GetElement(k - 1, j - 1) - 2 * res.GetElement(k - 1, j) + res.GetElement(k - 1, j + 1);
    };
    os << n << '\n';
    os << data.beta - data.alpha / data.h << ' ' << data.alpha / data.h << '\n';
    b.SetElement(0, 0, data.left_constraint(t_k));
    for (int j = 1; j < n - 1; ++j) {
        auto current_x = GetIthPoint(data.x_left, data.h, j);
        os << sigma * theta << ' ' << -(1 + 2 * sigma * theta) << ' ' << sigma * theta << '\n';
        b.SetElement(j, 0, -(res.GetElement(k - 1, j) + data.tau * data.f(current_x, t_k) +
                             (1 - theta) * sigma * getLowerRowSum(j)));
    }
    os << -data.gamma / data.h << ' ' << data.delta + data.gamma / data.h << '\n';
    b.SetElement(n - 1, 0, data.right_constraint(t_k));
    os >> a;
    return SweepMethod(a, b).GetColumn(0);
}

TMatrix CrankNicolsonSolution(const TaskData& data, double theta) {
    int number_of_rows = (data.t_top - data.t_bottom) / data.tau;
    int number_of_columns = (data.x_right - data.x_left) / data.h + 1;
    TMatrix result(number_of_rows, number_of_columns);
    for (int i = 0; i < number_of_columns; ++i) {
        result.SetElement(0, i, data.bottom_constraint(GetIthPoint(data.x_left, data.h, i)));
    }
    for (int i = 1; i < number_of_rows; ++i) {
        auto row = CrankNicolsonSolutionSingleStep(result, data, i, theta);
        for (int k = 0; k < number_of_columns; ++k) {
            result.SetElement(i, k, row[k]);
        }
    }
    return result;
}
