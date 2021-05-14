#include "../include/interpolators.h"

double lagrangeInterpolation(const std::vector<double>& i_xValues,
                             const std::vector<double>& i_yValues,
                             double i_point) {
    assert(i_xValues.size() == i_yValues.size());
    double result = 0;
    const double w = [&]() {
        double w = 1;
        for(double i_xValue : i_xValues) {
            w *= (i_point - i_xValue);
        }
        return w;
    }();
    auto getW = [&](int j) {
        double result = 1;
        for(int i = 0; i < i_xValues.size(); ++i) {
            if(i == j) {
                continue;
            }
            result *= (i_xValues[j] - i_xValues[i]);
        }
        return result;
    };
    for(int i = 0; i < i_xValues.size(); ++i) {
        result += i_yValues[i] * w / ((i_point - i_xValues[i]) * getW(i));
    }
    return result;
}

double newtonInterpolation(const std::vector<double> &i_xValues,
                           const std::vector<double> &i_yValues,
                           double i_point) {
    assert(i_xValues.size() == i_yValues.size());
    double result = 0;
    std::vector<std::vector<double>> differenceValues(i_xValues.size());
    for(double i_yValue : i_yValues) {
        differenceValues[0].push_back(i_yValue);
    }
    for(int i = 0; i < differenceValues.size() - 1; ++i) {
        for(int j = 0; j < differenceValues[i].size() - 1; ++j) {
            differenceValues[i + 1].push_back((differenceValues[i][j] - differenceValues[i][j + 1])
                                         / (i_xValues[j] - i_xValues[j + i + 1]));
        }
    }
    auto getProduct = [&](int j) {
        double result = 1;
        for(int i = 0; i < j; ++i) {
            result *= (i_point - i_xValues[i]);
        }
        return result;
    };
    for(int i = 0; i < differenceValues.size(); ++i) {
        result += getProduct(i) * differenceValues[i][0];
    }
    return result;
}

SplineInterpolator::SplineInterpolator(const std::vector<double> &i_xValues,
                                       const std::vector<double> &i_yValues) {
    assert(i_xValues.size() == i_yValues.size());
    const int n = i_yValues.size() - 1;
    TTridiagonalMatrix a;
    TMatrix b(n - 1, 1);
    std::stringstream os;
    auto getH = [&](int i) {
        return i_xValues[i] - i_xValues[i - 1];
    };
    os << n - 1 << '\n';
    os << 2 * (getH(1) + getH(2)) << ' ' << getH(2) << '\n';
    b.SetElement(0, 0, 3 * ((i_yValues[2] - i_yValues[1]) / getH(2) - (i_yValues[1] - i_yValues[0])) / getH(1));
    for(int i = 3; i <= n - 1; ++i) {
        os << getH(i - 1) << ' ' << 2 * (getH(i - 1) + getH(i)) << ' ' << getH(i) << '\n';
        b.SetElement(i - 2, 0, 3 * ((i_yValues[i] - i_yValues[i - 1]) / getH(i)
                                  - ((i_yValues[i - 1] - i_yValues[i - 2]) / getH(i - 1))));
    }
    os << getH(n - 1) << ' ' << 2 * (getH(n - 1) + getH(n)) << '\n';
    b.SetElement(n - 2, 0, 3 * ((i_yValues[n] - i_yValues[n - 1]) / getH(n)
                                             - (i_yValues[n - 1] - i_yValues[n - 2])) / getH(n - 1));
    os >> a;
    TMatrix c = SweepMethod(a, b);
    std::vector<double> cValues;
    cValues.push_back(0);
    for(int i = 0; i < n - 1; ++i) {
        cValues.push_back(c.GetElement(i, 0));
    }
    for(int i = 1; i <= n; ++i) {
        d_splines[i_xValues[i - 1]].push_back(i_yValues[i - 1]);
        if(i == n) {
            d_splines[i_xValues[i - 1]].push_back((i_yValues[i] - i_yValues[i - 1]) / getH(i) - 2. * getH(i) * cValues[i - 1] / 3.);
        } else {
            d_splines[i_xValues[i - 1]].push_back((i_yValues[i] - i_yValues[i - 1]) / getH(i) - getH(i)
                                                                                            * (2 * cValues[i - 1] + cValues[i]) / 3.);
        }
        d_splines[i_xValues[i - 1]].push_back(cValues[i - 1]);
        if(i == n) {
            d_splines[i_xValues[i - 1]].push_back(-cValues[i - 1] / (3. * getH(i)));
        } else {
            d_splines[i_xValues[i - 1]].push_back((cValues[i] - cValues[i - 1]) / (3 * getH(i)));
        }
    }
    for(int i = 0; i < 4; ++i) {
        d_splines[i_xValues[n]].push_back(0);
    }
}

double SplineInterpolator::getValue(double i_point) const {
    auto it = d_splines.lower_bound(i_point);
    if(it == d_splines.begin() || it == std::prev(d_splines.end())) {
        throw std::invalid_argument("Point is out of spline's range");
    }
    --it;
    double result = 0;
    auto getProduct = [&](double point, int power) {
        return pow((i_point - point), power);
    };
    for(int i = 0; i < 4; ++i) {
        result += getProduct(it->first, i) * it->second[i];
    }
    return result;
}

void SplineInterpolator::getSplineInfo(std::ostream& o_os) const {
    int count = 1;
    for(auto it = d_splines.begin();; ++it) {
        if(next(it) == d_splines.end()) {
            break;
        }
        o_os << "Spline #" << count << ", range: [" << it->first << ", " << std::next(it)->first << "], ";
        o_os << "a: " << it->second[0] << ", b: " << it->second[1] << ", c: " << it->second[2] << ", d:" << it->second[3];
        o_os << '\n';
        ++count;
    }
}

void SplineInterpolator::getSplineInfoWithoutText(std::ostream &o_os) const {
    for(auto it = d_splines.begin();; ++it) {
        if(next(it) == d_splines.end()) {
            break;
        }
        o_os << it->first << " " << std::next(it)->first << " ";
        o_os << it->second[0] << " " << it->second[1] << " " << it->second[2] << " " << it->second[3];
        o_os << '\n';
    }
}
