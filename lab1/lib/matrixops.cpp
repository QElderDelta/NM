#include "matrixops.h"

TMatrix BackwardSubstitution(const TMatrix &a, const TMatrix &b, size_t column_number) {
    TMatrix x(a.GetSize().number_of_rows, 1);
    double row_sum;
    for(int i = a.GetSize().number_of_rows - 1; i >= 0; --i) {
        row_sum = 0;
        for(size_t j = i + 1; j < a.GetSize().number_of_cols; ++j) {
            row_sum += a.GetElement(i, j) * x.GetElement(j, 0);
        }
        x.SetElement(i, 0, (b.GetElement(i, column_number) - row_sum) / a.GetElement(i, i));
    }
    return x;
}

TMatrix ForwardSubstitution(const TMatrix &a, const TMatrix &b, size_t column_number) {
    TMatrix x(a.GetSize().number_of_rows, 1);
    double row_sum;
    for(size_t i = 0; i < a.GetSize().number_of_rows; ++i) {
        row_sum = 0;
        for(size_t j = 0; j < i; ++j) {
            row_sum += a.GetElement(i, j) * x.GetElement(j, 0);
        }
        x.SetElement(i, 0, b.GetElement(i, column_number) - row_sum);
    }
    return x;
}

double GetDeterminantUsingLU(const TMatrix &l, const TMatrix &u, const TMatrix &p) {
    bool parity_is_odd = FindParityOfPermutation(p);
    double result = 1;
    double product = 1;
    for(size_t i = 0; i < l.GetSize().number_of_rows; ++i) {
        product *= l.GetElement(i, i);
    }
    result *= product;
    product = 1;
    for(size_t i = 0; i < l.GetSize().number_of_rows; ++i) {
        product *= u.GetElement(i, i);
    }
    result *= product;
    if(parity_is_odd) {
        result *= -1;
    }
    return result;
}

TMatrix InverseMatrixUsingLU(const TMatrix &l, const TMatrix &u, const TMatrix &p) {
    TMatrix z(l.GetSize());
    TMatrix e(l.GetSize());
    TMatrix x;
    e = p * e;
    for(size_t j = 0; j < l.GetSize().number_of_cols; ++j) {
        x = ForwardSubstitution(l, e, j);
        for(size_t i = 0; i < l.GetSize().number_of_rows; ++i) {
            z.SetElement(i, j, x.GetElement(i, 0));
        }
    }
    for(size_t j = 0; j < l.GetSize().number_of_cols; ++j) {
        x = BackwardSubstitution(u, z, j);
        for(size_t i = 0; i < l.GetSize().number_of_rows; ++i) {
            e.SetElement(i, j, x.GetElement(i, 0));
        }
    }
    return e;
}

TMatrix SolveLinearSystemUsingLU(const TMatrix &l, const TMatrix &u, const TMatrix &p, const TMatrix &b) {
    TMatrix z(l.GetSize().number_of_rows, 1);
    TMatrix x(l.GetSize().number_of_rows, 1);
    TMatrix b_permutation;
    b_permutation = p * b;
    z = ForwardSubstitution(l, b_permutation, 0);
    x = BackwardSubstitution(u, z, 0);
    return x;
}

bool FindParityOfPermutation(TMatrix p) {
    size_t permutation_count = 0;
    for(size_t i = 0; i < p.GetSize().number_of_rows; ++i) {
        if(p.GetElement(i, i) != 1) {
            ++permutation_count;
            p.SetElement(i, i, 1);
            for(size_t j = i + 1; j < p.GetSize().number_of_rows; ++j) {
                if(p.GetElement(j, i) == 1) {
                    p.SetElement(j, i, 0);
                }
            }
        }
    }
    return permutation_count % 2;
}

TMatrix SweepMethod(const TTridiagonalMatrix &a, const TMatrix &b) {
    size_t n = b.GetSize().number_of_rows;
    std::vector<std::vector<double>> coefs(n - 1);
    TMatrix x(n, 1);
    coefs[0].push_back(-a.GetElement(0, 1) / a.GetElement(0, 0));
    coefs[0].push_back(b.GetElement(0, 0) / a.GetElement(0, 0));
    for(size_t i = 1; i < n - 1; ++i) {
        coefs[i].push_back(-a.GetElement(i, 2) / (a.GetElement(i, 1) + a.GetElement(i, 0) * coefs[i - 1][0]));
        coefs[i].push_back((b.GetElement(i, 0) - a.GetElement(i, 0) * coefs[i - 1][1]) /
                            (a.GetElement(i, 1) + a.GetElement(i, 0) * coefs[i - 1][0]));
    }
    x.SetElement(n - 1, 0, (b.GetElement(n - 1, 0) - a.GetElement(n - 1, 0) * coefs[n - 1 - 1][1]) /
                           (a.GetElement(n - 1, 1) + a.GetElement(n - 1, 0) * coefs[n - 1 - 1][0]));
    for(int i = n - 2; i >= 0; --i) {
        x.SetElement(i, 0, coefs[i][0] * x.GetElement(i + 1, 0) + coefs[i][1]);
    }
    return x;
}

TMatrix SimpleIterationsMethod(const TMatrix &m, const TMatrix &b, double eps, std::ostream& log_stream) {
    TMatrix alpha(m.GetSize());
    TMatrix beta = b;
    for(size_t i = 0; i < m.GetSize().number_of_rows; ++i) {
        for(size_t j = 0; j < m.GetSize().number_of_cols; ++j) {
            if(i == j) {
                alpha.SetElement(i, j, 0);
            } else {
                alpha.SetElement(i, j, -m.GetElement(i, j) / m.GetElement(i, i));
            }
        }
        beta.SetElement(i, 0, b.GetElement(i, 0) / m.GetElement(i, i));
    }
    TMatrix x = beta;
    TMatrix x_prev;
    size_t iter_count = 0;
    double eps_k;
    log_stream << "Matrix alpha l1 norm: " << alpha.Getl1Norm() << '\n';
    bool use_alpha_norm = alpha.Getl1Norm() < 1;
    while(!iter_count || eps_k > eps) {
        ++iter_count;
        x_prev = x;
        x = beta + alpha * x_prev;
        if(use_alpha_norm) {
            eps_k = alpha.Getl1Norm() / (1 - alpha.Getl1Norm()) * (x - x_prev).Getl1Norm();
        } else {
            eps_k = (x - x_prev).Getl1Norm();
        }
    }
    log_stream << "Simple iterations method took " << iter_count << " iterations" << '\n';
    return x;
}

TMatrix SeidelMethod(const TMatrix &m, const TMatrix &b, double eps, std::ostream &log_stream) {
    TMatrix alpha(m.GetSize());
    TMatrix beta = b;
    for(size_t i = 0; i < m.GetSize().number_of_rows; ++i) {
        for(size_t j = 0; j < m.GetSize().number_of_cols; ++j) {
            if(i == j) {
                alpha.SetElement(i, j, 0);
            } else {
                alpha.SetElement(i, j, -m.GetElement(i, j) / m.GetElement(i, i));
            }
        }
        beta.SetElement(i, 0, b.GetElement(i, 0) / m.GetElement(i, i));
    }
    TMatrix c(m.GetSize());
    TMatrix d(m.GetSize());
    c.Clear();
    d.Clear();
    for(size_t i = 0; i < m.GetSize().number_of_rows; ++i) {
        for(size_t j = 0; j < m.GetSize().number_of_cols; ++j) {
            if(j >= i) {
                d.SetElement(i, j, alpha.GetElement(i, j));
            } else {
                c.SetElement(i, j, alpha.GetElement(i, j));
            }
        }
    }
    auto inverse = (TMatrix(m.GetSize()) - c).InverseMatrix();
    alpha = inverse * d;
    beta = inverse * beta;
    TMatrix x = beta;
    TMatrix x_prev;
    size_t iter_count = 0;
    double eps_k;
    log_stream << "Matrix alpha l1 norm: " << alpha.Getl1Norm() << '\n';
    bool use_alpha_norm = alpha.Getl1Norm() < 1;
    while(!iter_count || eps_k > eps) {
        ++iter_count;
        x_prev = x;
        x = beta + alpha * x_prev;
        if(use_alpha_norm) {
            eps_k = alpha.Getl1Norm() / (1 - alpha.Getl1Norm()) * (x - x_prev).Getl1Norm();
        } else {
            eps_k = (x - x_prev).Getl1Norm();
        }
    }
    log_stream << "Seidel's method took " << iter_count << " iterations" << '\n';
    return x;
}

std::pair<TMatrix, TMatrix> JacobiRotationMethod(const TMatrix &a, double eps, std::ostream& log_stream) {
    size_t iter_count = 0;
    TMatrix u(a.GetSize());
    TMatrix u_k(a.GetSize());
    TMatrix a_k = a;
    double phi_k;
    size_t n = a.GetSize().number_of_rows;
    size_t m = a.GetSize().number_of_cols;
    std::pair<size_t, size_t> max_element_pos;
    double max_element;
    while(std::sqrt(GetSumOfSquaredNonDiagonalElements(a_k)) > eps) {
        max_element = -1;
        max_element_pos = std::make_pair(0, 0);
        for(size_t i = 0; i < n; ++i) {
            for(size_t j = 0; j < m; ++j) {
                if(i != j && std::abs(a_k.GetElement(i, j)) > max_element) {
                    max_element = std::abs(a_k.GetElement(i, j));
                    max_element_pos = std::make_pair(i, j);
                }
            }
        }
        phi_k = 0.5 * std::atan(2 * a_k.GetElement(max_element_pos.first, max_element_pos.second) / (
                    a_k.GetElement(max_element_pos.first, max_element_pos.first)
                        - a_k.GetElement(max_element_pos.second, max_element_pos.second)
                ));
        u_k.SetElement(max_element_pos.first, max_element_pos.first, cos(phi_k));
        u_k.SetElement(max_element_pos.second, max_element_pos.second, cos(phi_k));
        u_k.SetElement(max_element_pos.first, max_element_pos.second, sin(phi_k));
        u_k.SetElement(max_element_pos.second, max_element_pos.first, -sin(phi_k));
        a_k = u_k * a_k;
        u_k.Transpose();
        a_k = a_k * u_k;
        u = u * u_k;
        u_k.SetIdentity();
        ++iter_count;
    }
    log_stream << "Jacobi's rotation method took " << iter_count << " iterations" << '\n';
    return std::make_pair(a_k, u);
}

double GetSumOfSquaredNonDiagonalElements(const TMatrix &a) {
    double res = 0;
    for(size_t i = 0; i < a.GetSize().number_of_rows; ++i) {
        for(size_t j = 0; j < a.GetSize().number_of_cols; ++j) {
            if(i != j) {
                res += a.GetElement(i, j) * a.GetElement(i, j);
            }
        }
    }
    return res;
}
