#pragma once

#include "function.h"

#include <iostream>
#include <cassert>

double newtonsMethod(ISingleArgumentFunction& i_function, double i_x0, double i_eps,
                     std::ostream& log_stream = std::cout);

double simpleIterationsMethod(ISingleArgumentFunction& i_function, double i_x0, double i_q, double i_eps,
                              std::ostream& log_stream = std::cout);