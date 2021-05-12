#include "../include/methods.h"

double newtonsMethod(ISingleArgumentFunction& i_function, double i_x0, double i_eps, std::ostream& log_stream) {
    int numberOfIterations = 1;
    double xk = i_x0 - i_function.getValue(i_x0) / i_function.firstDerivative(i_x0);
    while(std::abs(xk - i_x0) > i_eps) {
        ++numberOfIterations;
        i_x0 = xk;
        xk = i_x0 - i_function.getValue(i_x0) / i_function.firstDerivative(i_x0);
    }
    log_stream << "Newton's method took " << numberOfIterations << " iterations" << '\n';
    return xk;
}

double simpleIterationsMethod(ISingleArgumentFunction& i_function, double i_x0, double i_q, double i_eps,
                              std::ostream& log_stream) {
    assert(i_q < 1);
    int numberOfIterations = 1;
    double xk = i_function.getValue(i_x0);
    while(i_q * std::abs(xk - i_x0) / (1 - i_q) > i_eps) {
        ++numberOfIterations;
        i_x0 = xk;
        xk = i_function.getValue(i_x0);
    }
    log_stream << "Simple iterations method took " << numberOfIterations << " iterations" << '\n';
    return xk;
}