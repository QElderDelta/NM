#include <2d_eq_solver.h>

#include <sstream>

namespace {
    double getIthPoint(double start, double h, int i) {
        return start + h * i;
    }

    void fillGrid(std::vector<TMatrix>& grid, const TaskData& data) {
        int k = grid.size();
        int n = grid[0].GetSize().number_of_rows;
        int m = grid[0].GetSize().number_of_cols;
        double curr_t, curr_x, curr_y;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                grid[0].SetElement(i, j, data.first_layer_constraint(getIthPoint(data.x_left, data.hx, j),
                                                                     getIthPoint(data.y_bottom, data.hy, i)));
            }
        }
        for (int t = 1; t < k; ++t) {
            curr_t = getIthPoint(data.t_bottom, data.tau, t);
            for (int i = 0; i < n; ++i) {
                curr_y = getIthPoint(data.y_bottom, data.hy, i);
                grid[t].SetElement(i, 0, data.left_constraint(curr_y, curr_t));
                grid[t].SetElement(i, m - 1, data.right_constraint(curr_y, curr_t));
            }
            for (int j = 0; j < n; ++j) {
                curr_x = getIthPoint(data.x_left, data.hx, j);
                grid[t].SetElement(0, j, data.bottom_constraint(curr_x, curr_t));
                grid[t].SetElement(n - 1, j, data.top_constraint(curr_x, curr_t));
            }
        }
    }
}

std::vector<TMatrix> varyingDirectionsMethod(const TaskData& data) {
    int number_of_rows = (data.y_top - data.y_bottom) / data.hy + 1;
    int number_of_columns = (data.x_right - data.x_left) / data.hx + 1;
    int number_of_layers = (data.t_top - data.t_bottom) / data.tau + 1;
    std::vector<TMatrix> result(number_of_layers, TMatrix(number_of_rows, number_of_columns));
    std::stringstream os;
    TMatrix temp_lvl;
    TTridiagonalMatrix a;
    TMatrix b(number_of_columns, 1), c;
    double curr_t, curr_y, curr_x;
    double hx_2 = pow(data.hx, 2);
    double hy_2 = pow(data.hy, 2);
    double u_sum;
    fillGrid(result, data);
    for (int k = 1; k < number_of_layers; ++k) {
        curr_t = getIthPoint(data.t_bottom, data.tau, k);
        temp_lvl = TMatrix(number_of_rows, number_of_columns);
        for (int i = 1; i < number_of_rows - 1; ++i) {
            curr_y = getIthPoint(data.y_bottom, data.hy, i);
            os.clear();
            os << number_of_columns << '\n';
            os << "1 0" << '\n';
            b.SetElement(0, 0, data.left_constraint(curr_y, curr_t - data.tau / 2.));
            for (int j = 1; j < number_of_columns - 1; ++j) {
                curr_x = getIthPoint(data.x_left, data.hx, j);
                u_sum = 0;
                os << -data.a * data.tau / (2 * hx_2) << ' ';
                os << (1 + data.a * data.tau / hx_2) << ' ';
                os << -data.a * data.tau / (2 * hx_2) << ' ';
                u_sum += result[k - 1].GetElement(i - 1, j);
                u_sum += -2 * result[k - 1].GetElement(i, j);
                u_sum += result[k - 1].GetElement(i + 1, j);
                b.SetElement(j, 0, result[k - 1].GetElement(i, j) +
                                   data.tau *
                                   (data.b * u_sum / hy_2 + data.left_part(curr_x, curr_y, curr_t - data.tau / 2.)) /
                                   2.);
            }
            os << "0 1" << '\n';
            b.SetElement(number_of_columns - 1, 0, data.right_constraint(curr_y, curr_t - data.tau / 2.));
            os >> a;
            c = SweepMethod(a, b);
            for (int j = 0; j < c.GetSize().number_of_rows; ++j) {
                temp_lvl.SetElement(i, j, c.GetElement(j, 0));
            }
        }

        b = TMatrix(number_of_rows, 1);
        for (int j = 1; j < number_of_columns - 1; ++j) {
            curr_x = getIthPoint(data.x_left, data.hx, j);
            os.clear();
            os << number_of_rows << '\n';
            os << "1 0" << '\n';
            b.SetElement(0, 0, data.bottom_constraint(curr_x, curr_t));
            for (int i = 1; i < number_of_rows - 1; ++i) {
                u_sum = 0;
                os << -data.b * data.tau / (2 * hy_2) << ' ';
                os << 1 + data.b * data.tau / hy_2 << ' ';
                os << -data.b * data.tau / (2 * hy_2) << ' ';
                u_sum += temp_lvl.GetElement(i, j - 1);
                u_sum += -2. * temp_lvl.GetElement(i, j);
                u_sum += temp_lvl.GetElement(i, j + 1);
                b.SetElement(i, 0, temp_lvl.GetElement(i, j) + data.tau * (data.a * u_sum / hx_2 +
                                                                           data.left_part(curr_x, curr_y,
                                                                                          curr_t - data.tau / 2.)) /
                                                               2.);
            }
            os << "0 1" << '\n';
            b.SetElement(number_of_rows - 1, 0, data.top_constraint(curr_x, curr_t));
            os >> a;
            c = SweepMethod(a, b);
            for (int i = 0; i < c.GetSize().number_of_rows; ++i) {
                result[k].SetElement(i, j, c.GetElement(i, 0));
            }
        }
    }
    return result;
}

std::vector<TMatrix> fractionalStepsMethod(const TaskData& data) {
    int number_of_rows = (data.y_top - data.y_bottom) / data.hy + 1;
    int number_of_columns = (data.x_right - data.x_left) / data.hx + 1;
    int number_of_layers = (data.t_top - data.t_bottom) / data.tau + 1;
    std::vector<TMatrix> result(number_of_layers, TMatrix(number_of_rows, number_of_columns));
    std::stringstream os;
    TMatrix temp_lvl;
    TTridiagonalMatrix a;
    TMatrix b(number_of_columns, 1), c;
    double curr_t, curr_y, curr_x;
    double hx_2 = pow(data.hx, 2);
    double hy_2 = pow(data.hy, 2);
    fillGrid(result, data);
    for (int k = 1; k < number_of_layers; ++k) {
        curr_t = getIthPoint(data.t_bottom, data.tau, k);
        temp_lvl = TMatrix(number_of_rows, number_of_columns);
        for (int i = 1; i < number_of_rows - 1; ++i) {
            curr_y = getIthPoint(data.y_bottom, data.hy, i);
            os.clear();
            os << number_of_columns << '\n';
            os << "1 0" << '\n';
            b.SetElement(0, 0, data.left_constraint(curr_y, curr_t - data.tau / 2.));
            for (int j = 1; j < number_of_columns - 1; ++j) {
                curr_x = getIthPoint(data.x_left, data.hx, j);
                os << -data.a * data.tau / hx_2 << ' ';
                os << (1 + 2 * data.a * data.tau / hx_2) << ' ';
                os << -data.a * data.tau / hx_2 << ' ';
                b.SetElement(j, 0, result[k - 1].GetElement(i, j) +
                                   data.tau * data.left_part(curr_x, curr_y, curr_t - data.tau) / 2.);
            }
            os << "0 1" << '\n';
            b.SetElement(number_of_columns - 1, 0, data.right_constraint(curr_y, curr_t - data.tau / 2.));
            os >> a;
            c = SweepMethod(a, b);
            for (int j = 0; j < c.GetSize().number_of_rows; ++j) {
                temp_lvl.SetElement(i, j, c.GetElement(j, 0));
            }
        }

        b = TMatrix(number_of_rows, 1);
        for (int j = 1; j < number_of_columns - 1; ++j) {
            curr_x = getIthPoint(data.x_left, data.hx, j);
            os.clear();
            os << number_of_rows << '\n';
            os << "1 0" << '\n';
            b.SetElement(0, 0, data.bottom_constraint(curr_x, curr_t));
            for (int i = 1; i < number_of_rows - 1; ++i) {
                os << -data.b * data.tau / hy_2 << ' ';
                os << 1 + 2 * data.b * data.tau / hy_2 << ' ';
                os << -data.b * data.tau / hy_2 << ' ';
                b.SetElement(i, 0, temp_lvl.GetElement(i, j) + data.tau * data.left_part(curr_x, curr_y, curr_t) / 2.);
            }
            os << "0 1" << '\n';
            b.SetElement(number_of_rows - 1, 0, data.top_constraint(curr_x, curr_t));
            os >> a;
            c = SweepMethod(a, b);
            for (int i = 0; i < c.GetSize().number_of_rows; ++i) {
                result[k].SetElement(i, j, c.GetElement(i, 0));
            }
        }
    }
    return result;
}