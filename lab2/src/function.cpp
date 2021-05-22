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

double SecondTaskFunction1::getValue(double i_x1, double i_x2) {
    return i_x1 - cos(i_x2) - 1;
}

double SecondTaskFunction1::getFirstDerivativeByFirstArgument(double i_x1, double i_x2) {
    return 1;
}

double SecondTaskFunction1::getFirstDerivativeBySecondArgument(double i_x1, double i_x2) {
    return sin(i_x2);
}

double SecondTaskFunction1::getValue(std::pair<double, double> i_point) {
    return this->getValue(i_point.first, i_point.second);
}

double SecondTaskFunction2::getValue(double i_x1, double i_x2) {
    return i_x2 - log10(i_x1 + 1) - 2;
}

double SecondTaskFunction2::getFirstDerivativeByFirstArgument(double i_x1, double i_x2) {
    return -log(i_x1) / log(10);
}

double SecondTaskFunction2::getFirstDerivativeBySecondArgument(double i_x1, double i_x2) {
    return 1;
}

double SecondTaskFunction2::getValue(std::pair<double, double> i_point) {
    return this->getValue(i_point.first, i_point.second);
}

double SecondTaskFunction1Transformed::getValue(double i_x1, double i_x2) {
    return 1 + cos(i_x2);
}

double SecondTaskFunction2Transformed::getValue(double i_x1, double i_x2) {
    return 2 + log10(i_x1 + 1);
}
