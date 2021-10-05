#pragma once

#include <matrix.h>

#include <functional>
#include <vector>

using ConstraintFunc = std::function<double(double)>;

using FreeFunc = std::function<double(double, double)>;

struct TaskData {
    ConstraintFunc left_constraint;
    double alpha, beta;
    ConstraintFunc right_constraint;
    double gamma, delta;
    ConstraintFunc bottom_constraint;
    double tau, h, a;
    double x_left, x_right;
    double t_bottom, t_top;
    FreeFunc f;
};

double GetIthPoint(double start, double h, int i);

std::vector<double>
ExplicitSchemaSolutionSingleStep(const TMatrix& res, const TaskData& data, int k, std::optional<double> theta = {});

TMatrix ExplicitSchemaSolution(const TaskData& data);

std::vector<double> ImplicitSchemaSolutionSingleStep(const TMatrix& res, const TaskData& data, int k,
                                                     std::optional<std::vector<double>> to_add = {},
                                                     std::optional<double> theta = {});

TMatrix ImplicitSchemaSolution(const TaskData& data);

std::vector<double> CrankNicolsonSolutionSingleStep(const TMatrix& res, const TaskData& data, int k, double theta);

TMatrix CrankNicolsonSolution(const TaskData& data, double theta = 0.5);