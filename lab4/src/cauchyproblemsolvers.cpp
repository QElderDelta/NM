#include "../include/cauchyproblemsolvers.h"

std::vector<std::pair<double, double>> eulerMethod(ITwoArgumentFunction &i_f,
                                                   double i_a,
                                                   double i_b,
                                                   double i_h,
                                                   double i_leftConstraint) {
    std::vector<std::pair<double, double>> result;
    while(i_a <= i_b) {
        result.emplace_back(i_a, i_leftConstraint);
        i_leftConstraint = i_leftConstraint + i_h * i_f.getValue(i_a, i_leftConstraint);
        i_a += i_h;
    }
    return result;
}

std::vector<std::pair<double, double>> rungeKuttaFourthOrderMethod(ITwoArgumentFunction &i_f,
                                                                   double i_a,
                                                                   double i_b,
                                                                   double i_h,
                                                                   double i_leftConstraint) {
    std::vector<std::pair<double, double>> result;
    double k1, k2, k3, k4;
    while(i_a <= i_b) {
        result.emplace_back(i_a, i_leftConstraint);
        k1 = i_h * i_f.getValue(i_a, i_leftConstraint);
        k2 = i_h * i_f.getValue(i_a + 0.5 * i_h, i_leftConstraint + 0.5 * k1);
        k3 = i_h * i_f.getValue(i_a + 0.5 * i_h, i_leftConstraint + 0.5 * k2);
        k4 = i_h * i_f.getValue(i_a + i_h, i_leftConstraint + k3);
        i_leftConstraint = i_leftConstraint + (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
        i_a += i_h;
    }
    return result;
}

std::vector<std::pair<double, double>> adamsFourthOrderMethod(ITwoArgumentFunction &i_f,
                                                              double i_a,
                                                              double i_b,
                                                              double i_h,
                                                              double i_leftConstraint) {
    std::vector<std::pair<double, double>> result;
    auto first_steps = rungeKuttaFourthOrderMethod(i_f, i_a, std::min(i_a +  3 * i_h, i_b), i_h, i_leftConstraint);
    result.insert(result.end(), std::make_move_iterator(first_steps.begin()),
                                std::make_move_iterator(first_steps.end()));
    i_a = i_a + 4 * i_h;
    auto getSum = [&i_f, &result](int counter) {
        return (55 * i_f.getValue(result[counter]) - 59 * i_f.getValue(result[counter - 1]) +
                37 * i_f.getValue(result[counter - 2]) - 9 * i_f.getValue(result[counter - 3]));
    };
    int counter = result.size() - 1;
    i_leftConstraint = result.back().second + i_h * getSum(counter) / 24.;
    while(i_a <= i_b) {
        ++counter;
        result.emplace_back(i_a, i_leftConstraint);
        i_leftConstraint = i_leftConstraint + i_h * getSum(counter) / 24.;
        i_a += i_h;
    }
    return result;
}

std::vector<Point> solveSecondOrderUsingEulerMethod(Function &i_f,
                                                    double i_a,
                                                    double i_b,
                                                    double i_h,
                                                    double i_yLeftConstraint,
                                                    double i_zLeftConstraint) {
    std::vector<Point> result;
    while(i_a <= i_b) {
        result.emplace_back(i_a, i_yLeftConstraint, i_zLeftConstraint);
        i_yLeftConstraint = i_yLeftConstraint + i_h * i_zLeftConstraint;
        i_zLeftConstraint = i_zLeftConstraint + i_h * i_f.getValue({i_a, i_yLeftConstraint, i_zLeftConstraint});
        i_a += i_h;
    }
    return result;
}

std::vector<Point> solveSecondOrderUsingRungeKuttaFourthOrderMethod(Function &i_f,
                                                                    double i_a,
                                                                    double i_b,
                                                                    double i_h,
                                                                    double i_yLeftConstraint,
                                                                    double i_zLeftConstraint) {
    std::vector<Point> result;
    double k1, k2, k3, k4, l1, l2, l3, l4;
    while(i_a <= i_b) {
        result.emplace_back(i_a, i_yLeftConstraint, i_zLeftConstraint);
        k1 = i_h * i_zLeftConstraint;
        l1 = i_h * i_f.getValue({i_a, i_yLeftConstraint, i_zLeftConstraint});
        k2 = i_h * (i_zLeftConstraint + 0.5 * l1);
        l2 = i_h * i_f.getValue({i_a + 0.5 * i_h, i_yLeftConstraint + 0.5 * k1, i_zLeftConstraint + 0.5 * l1});
        k3 = i_h * (i_zLeftConstraint + 0.5 * l2);
        l3 = i_h * i_f.getValue({i_a + 0.5 * i_h, i_yLeftConstraint + 0.5 * k2, i_zLeftConstraint + 0.5 * l2});
        k4 = i_h * (i_zLeftConstraint + l3);
        l4 = i_h * i_f.getValue({i_a + i_h, i_yLeftConstraint + k3, i_zLeftConstraint + l3});
        i_yLeftConstraint = i_yLeftConstraint + (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
        i_zLeftConstraint = i_zLeftConstraint + (l1 + 2 * l2 + 2 * l3 + l4) / 6.;
        i_a += i_h;
    }
    return result;
}

std::vector<Point> solveSecondOrderUsingAdamsFourthOrderMethod(Function &i_f,
                                                               double i_a,
                                                               double i_b,
                                                               double i_h,
                                                               double i_yLeftConstraint,
                                                               double i_zLeftConstraint) {
    std::vector<Point> result;
    auto first_steps = solveSecondOrderUsingRungeKuttaFourthOrderMethod(i_f,
                                                                        i_a,
                                                                        std::min(i_a +  3 * i_h, i_b),
                                                                        i_h,
                                                                        i_yLeftConstraint,
                                                                        i_zLeftConstraint);
    result.insert(result.end(), std::make_move_iterator(first_steps.begin()),
                  std::make_move_iterator(first_steps.end()));
    i_a = i_a + 4 * i_h;
    auto makeVector = [&result](int counter) {
        return std::vector<double>{result[counter].x, result[counter].y, result[counter].z};
    };
    auto getYSum = [&result](int counter) {
        return (55 * result[counter].z - 59 * result[counter - 1].z +
                37 * result[counter - 2].z - 9 * result[counter - 3].z);
    };
    auto getZSum = [&i_f, &makeVector](int counter) {
        return (55 * i_f.getValue(makeVector(counter)) - 59 * i_f.getValue(makeVector(counter - 1)) +
                37 * i_f.getValue(makeVector(counter - 2)) - 9 * i_f.getValue(makeVector(counter - 3)));
    };
    int counter = result.size() - 1;
    i_yLeftConstraint = result.back().y + i_h * getYSum(counter) / 24.;
    i_zLeftConstraint = result.back().z + i_h * getZSum(counter) / 24.;
    while(i_a <= i_b) {
        ++counter;
        result.emplace_back(i_a, i_yLeftConstraint, i_zLeftConstraint);
        i_yLeftConstraint = result.back().y + i_h * getYSum(counter) / 24.;
        i_zLeftConstraint = result.back().z + i_h * getZSum(counter) / 24.;
        i_a += i_h;
    }
    return result;
}


Function::Function(int order) : d_order(order) {}

int Function::getOrder() const {
    return d_order;
}
