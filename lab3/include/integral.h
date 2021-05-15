#pragma once

#include <cassert>
#include <cmath>
#include <functional>

double integrateUsingRectangleMethod(const std::function<double(double)>& i_f,
                                     double i_start,
                                     double i_end,
                                     double i_step);

double integrateUsingTrapezoidMethod(const std::function<double(double)>& i_f,
                                     double i_start,
                                     double i_end,
                                     double i_step);

double integrateUsingSimpsonMethod(const std::function<double(double)>& i_f,
                                   double i_start,
                                   double i_end,
                                   double i_step);

double rungeRombergMethod(double i_fh, double i_fkh, double k, double p);