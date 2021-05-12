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

std::pair<double, double>
newtonsMethod(ITwoArgumentFunction &i_f1, ITwoArgumentFunction &i_f2, std::pair<double, double> i_x0,
              double i_eps, std::ostream &log_stream) {
    int numberOfIterations = 1;
    TMatrix x0(2, 1);
    x0.SetElement(0, 0, i_x0.first);
    x0.SetElement(1, 0, i_x0.second);
    TMatrix f(2, 1);
    f.SetElement(0, 0, i_f1.getValue(x0.GetElement(0, 0), x0.GetElement(1, 0)));
    f.SetElement(1, 0, i_f2.getValue(x0.GetElement(0, 0), x0.GetElement(1, 0)));
    TMatrix J(2, 2);
    J.SetElement(0, 0, i_f1.getFirstDerivativeByFirstArgument(x0.GetElement(0, 0),
                                                              x0.GetElement(1, 0)));
    J.SetElement(0, 1, i_f1.getFirstDerivativeBySecondArgument(x0.GetElement(0, 0),
                                                               x0.GetElement(1, 0)));
    J.SetElement(1, 0, i_f2.getFirstDerivativeByFirstArgument(x0.GetElement(0, 0),
                                                              x0.GetElement(1, 0)));
    J.SetElement(1, 1, i_f2.getFirstDerivativeBySecondArgument(x0.GetElement(0, 0),
                                                               x0.GetElement(1, 0)));
    TMatrix xk = x0 - J.InverseMatrix() * f;
    while((xk - x0).Getl1Norm() > i_eps) {
        ++numberOfIterations;
        x0 = xk;
        f.SetElement(0, 0, i_f1.getValue(x0.GetElement(0, 0), x0.GetElement(1, 0)));
        f.SetElement(1, 0, i_f2.getValue(x0.GetElement(0, 0), x0.GetElement(1, 0)));
        J.SetElement(0, 0, i_f1.getFirstDerivativeByFirstArgument(x0.GetElement(0, 0),
                                                                  x0.GetElement(1, 0)));
        J.SetElement(0, 1, i_f1.getFirstDerivativeBySecondArgument(x0.GetElement(0, 0),
                                                                   x0.GetElement(1, 0)));
        J.SetElement(1, 0, i_f2.getFirstDerivativeByFirstArgument(x0.GetElement(0, 0),
                                                                  x0.GetElement(1, 0)));
        J.SetElement(1, 1, i_f2.getFirstDerivativeBySecondArgument(x0.GetElement(0, 0),
                                                                   x0.GetElement(1, 0)));
        xk = x0 - J.InverseMatrix() * f;
    }
    log_stream << "Newton's method took " << numberOfIterations << " iterations" << '\n';
    return {xk.GetElement(0, 0), xk.GetElement(1, 0)};
}

std::pair<double, double>
simpleIterationsMethod(ITwoArgumentFunction &i_f1, ITwoArgumentFunction &i_f2, std::pair<double, double> i_x0,
                       double i_q, double i_eps, std::ostream &log_stream) {
    assert(i_q < 1);
    int numberOfIterations = 1;
    TMatrix x0(2, 1);
    x0.SetElement(0, 0, i_x0.first);
    x0.SetElement(1, 0, i_x0.second);
    TMatrix xk(2, 1);
    xk.SetElement(0, 0, i_f1.getValue(x0.GetElement(0, 0), x0.GetElement(1, 0)));
    xk.SetElement(1, 0, i_f2.getValue(x0.GetElement(0, 0), x0.GetElement(1, 0)));
    while(i_q * (xk - x0).Getl1Norm() / (1 - i_q) > i_eps) {
        ++numberOfIterations;
        x0 = xk;
        xk.SetElement(0, 0, i_f1.getValue(x0.GetElement(0, 0), x0.GetElement(1, 0)));
        xk.SetElement(1, 0, i_f2.getValue(x0.GetElement(0, 0), x0.GetElement(1, 0)));
    }
    log_stream << "Simple iterations method took " << numberOfIterations << " iterations" << '\n';
    return {xk.GetElement(0, 0), xk.GetElement(1, 0)};
}
