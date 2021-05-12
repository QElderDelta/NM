#pragma once

#include "../../lab1/lib/matrix.h"
#include "function.h"

#include <iostream>
#include <cassert>

double newtonsMethod(ISingleArgumentFunction& i_function, double i_x0, double i_eps,
                     std::ostream& log_stream = std::cout);

double simpleIterationsMethod(ISingleArgumentFunction& i_function, double i_x0, double i_q, double i_eps,
                              std::ostream& log_stream = std::cout);

std::pair<double, double> newtonsMethod(ITwoArgumentFunction& i_f1, ITwoArgumentFunction& i_f2,
                                        std::pair<double, double> i_x0, double i_eps,
                                        std::ostream& log_stream = std::cout);

std::pair<double, double> simpleIterationsMethod(ITwoArgumentFunction& i_f1, ITwoArgumentFunction& i_f2,
                                                 std::pair<double, double> i_x0, double i_q, double i_eps,
                                                 std::ostream& log_stream = std::cout);