#pragma once

#include <matrix.h>
#include <matrixops.h>

#include <functional>

using ConstraintFunc = std::function<double(double)>;

using FreeFunc = std::function<double(double, double)>;

//d^2u/dt^2 + e * du/dt = a^2 * d^2u/dx^2 + b * du/dx + c * u + f
//alpha * du(0,t)/dx + beta * u(0,t) = left_constraint
//gamma * du(l,t)/dx + delta * u(l,t) = right_constraint
//u(x,0) = bottom_constraint
//du(x,0)/dt = bottom_constraint_derivative
struct TaskData {
    ConstraintFunc left_constraint;
    double alpha, beta;
    ConstraintFunc right_constraint;
    double gamma, delta;
    ConstraintFunc bottom_constraint;
    ConstraintFunc bottom_constraint_derivative;
    ConstraintFunc bottom_constraint_second_derivative;
    double tau, h;
    double a, b, c, e;
    double x_left, x_right;
    double t_bottom, t_top;
    FreeFunc f;
};

TMatrix solveHyperbolicEquation(const TaskData& data, double s);