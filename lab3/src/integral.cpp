#include "../include/integral.h"

double integrateUsingRectangleMethod(const std::function<double(double)> &i_f,
                                     double i_start,
                                     double i_end,
                                     double i_step) {
    double result = 0;
    double current = i_start;
    while(current + i_step <= i_end) {
        result += i_f(current + i_step / 2.);
        current += i_step;
    }
    return i_step * result;
}

double integrateUsingTrapezoidMethod(const std::function<double(double)> &i_f,
                                     double i_start,
                                     double i_end,
                                     double i_step) {
    double result = 0;
    double current = i_start;
    result += i_f(current) / 2.;
    current += i_step;
    while(true) {
        if(current + i_step > i_end) {
            result += i_f(current) / 2.;
            break;
        } else {
            result += i_f(current);
            current += i_step;
        }
    }
    return i_step * result;
}

double integrateUsingSimpsonMethod(const std::function<double(double)> &i_f,
                                   double i_start,
                                   double i_end,
                                   double i_step) {
    double result = 0;
    double current = i_start;
    bool multiplyBy4 = true;
    result += i_f(current);
    current += i_step;
    while(true) {
        if(current + i_step > i_end) {
            result += i_f(current);
            break;
        } else {
            if(multiplyBy4) {
                result += 4. * i_f(current);
            } else {
                result += 2. * i_f(current);
            }
            current += i_step;
            multiplyBy4 = !multiplyBy4;
        }
    }
    return i_step * result / 3.;
}

double rungeRombergMethod(double i_fh, double i_fkh, double k, double p) {
    assert(k != 1 && k > 0);
    return i_fh + (i_fh - i_fkh) / (pow(k, p) - 1);
}
