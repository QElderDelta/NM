#include <elliptic_eq_solver.h>

namespace {
    double getIthPoint(double start, double h, int i) {
        return start + h * i;
    }

    TMatrix fillGrid(const TaskData& data) {
        double x_range = data.x_right - data.x_left;
        int number_of_rows = (data.y_top - data.y_bottom) / data.hy;
        int number_of_columns = x_range / data.hx + 1;
        TMatrix result(number_of_rows, number_of_columns);
        double curr_point, u_left, u_right;
        for (int i = 0; i < number_of_rows; ++i) {
            curr_point = getIthPoint(data.y_bottom, data.hy, i);
            result.SetElement(i, 0, data.left_constraint(curr_point));
            result.SetElement(i, number_of_columns - 1, data.right_constraint(curr_point));
        }
        for (int i = 0; i < number_of_columns; ++i) {
            curr_point = getIthPoint(data.x_left, data.hx, i);
            result.SetElement(0, i, data.bottom_constraint(curr_point));
            result.SetElement(number_of_rows - 1, i, data.top_constraint(curr_point));
        }
        for (int i = 1; i < number_of_rows - 1; ++i) {
            u_left = result.GetElement(i, 0);
            u_right = result.GetElement(i, number_of_columns - 1);
            for (int j = 1; j < number_of_columns - 1; ++j) {
                result.SetElement(i, j, u_left +
                                        (u_right - u_left) * (getIthPoint(data.x_left, data.hx, j) - data.x_left) /
                                        x_range);
            }
        }

        return result;
    }
}

TMatrix liebmannMethod(const TaskData& data, std::ostream* log_stream) {
    auto prev = fillGrid(data);
    double max_diff = -1;
    double u_curr;
    int number_of_rows = prev.GetSize().number_of_rows;
    int number_of_cols = prev.GetSize().number_of_cols;
    TMatrix curr = prev;
    double hx_2 = pow(data.hx, 2);
    double hy_2 = pow(data.hy, 2);
    while (max_diff == -1 || max_diff > data.eps) {
        if (max_diff != -1) {
            prev = curr;
            if (log_stream) {
                *log_stream << max_diff << ' ';
            }
        }
        max_diff = -1;
        for (int i = 1; i < number_of_rows - 1; ++i) {
            for (int j = 1; j < number_of_cols - 1; ++j) {
                u_curr = 0;
                u_curr += hx_2 * (prev.GetElement(i - 1, j) + prev.GetElement(i + 1, j));
                u_curr += hy_2 * (prev.GetElement(i, j - 1) + prev.GetElement(i, j + 1));
                u_curr /= (2 * hx_2 + 2 * hy_2 - 2 * hx_2 * hy_2);
                curr.SetElement(i, j, u_curr);
                if (std::fabs(prev.GetElement(i, j) - u_curr) > max_diff) {
                    max_diff = std::fabs(prev.GetElement(i, j) - u_curr);
                }
            }
        }
    }
    return curr;
}

TMatrix relaxationMethod(const TaskData& data, double tau, std::ostream* log_stream) {
    auto prev = fillGrid(data);
    double max_diff = -1;
    double u_curr;
    int number_of_rows = prev.GetSize().number_of_rows;
    int number_of_cols = prev.GetSize().number_of_cols;
    TMatrix curr = prev;
    double hx_2 = pow(data.hx, 2);
    double hy_2 = pow(data.hy, 2);
    while (max_diff == -1 || max_diff > data.eps) {
        if (max_diff != -1) {
            prev = curr;
            if (log_stream) {
                *log_stream << max_diff << ' ';
            }
        }
        max_diff = -1;
        for (int i = 1; i < number_of_rows - 1; ++i) {
            for (int j = 1; j < number_of_cols - 1; ++j) {
                u_curr = 0;
                u_curr += hx_2 * (prev.GetElement(i - 1, j) + prev.GetElement(i + 1, j));
                u_curr += hy_2 * (prev.GetElement(i, j - 1) + prev.GetElement(i, j + 1));
                u_curr -= hx_2 * hy_2 * data.left_eq_part(prev.GetElement(i, j));
                u_curr /= 2 * (hx_2 + hy_2);
                u_curr -= (1 - 1 / tau) * prev.GetElement(i, j);
                u_curr *= tau;
                curr.SetElement(i, j, u_curr);
                if (std::fabs(prev.GetElement(i, j) - u_curr) > max_diff) {
                    max_diff = std::fabs(prev.GetElement(i, j) - u_curr);
                }
            }
        }
    }
    return curr;
}

TMatrix seidelMethod(const TaskData& data, std::ostream* log_stream) {
    return relaxationMethod(data, 1, log_stream);
}