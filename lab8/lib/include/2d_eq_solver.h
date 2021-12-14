#pragma once

#include <matrixops.h>

#include <functional>
#include <optional>

using ConstraintFunc = std::function<double(double, double)>;
using FreeFunc = std::function<double(double, double, double)>;

struct TaskData {
    double hx, hy, tau;
    double x_left, x_right;
    double y_bottom, y_top;
    double t_top, t_bottom;
    double a, b;
    ConstraintFunc left_constraint, right_constraint;
    ConstraintFunc bottom_constraint, top_constraint;
    ConstraintFunc first_layer_constraint;
    FreeFunc left_part;
};

std::vector<TMatrix> varyingDirectionsMethod(const TaskData& data);

std::vector<TMatrix> fractionalStepsMethod(const TaskData& data);
