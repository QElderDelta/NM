#include "../include/function.h"

double TaskFunction::getValue(double i_x) {
    return cos(i_x) + 0.25 * i_x - 0.5;
}

double TaskFunction::firstDerivative(double i_x) {
    return -sin(i_x) + 0.25;
}

double TransformedFunction::getValue(double i_x) {
    return acos(-0.25 * i_x + 0.5);
}
