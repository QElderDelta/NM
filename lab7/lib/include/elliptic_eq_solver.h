#pragma once

#include <matrixops.h>

#include <functional>
#include <optional>

using ConstraintFunc = std::function<double(double)>;

struct TaskData {
    double hx, hy;
    double x_left, x_right;
    double y_bottom, y_top;
    double eps;
    ConstraintFunc left_constraint, right_constraint;
    ConstraintFunc bottom_constraint, top_constraint;
    ConstraintFunc left_eq_part;
};

TMatrix liebmannMethod(const TaskData& data, std::ostream* log_stream = nullptr);

TMatrix relaxationMethod(const TaskData& data, double tau, std::ostream* log_stream = nullptr);

TMatrix seidelMethod(const TaskData& data, std::ostream* log_stream = nullptr);
