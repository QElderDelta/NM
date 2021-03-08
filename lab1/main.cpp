#include <fstream>
#include <iostream>
#include "lib/matrix.h"

int main() {
    TMatrix m1, m2;
    std::ifstream is("../test.txt");
    is >> m1 >> m2;
    auto [l, u, p] = m1.LUDecomposition();
    //std::cout << l << '\n' << u << '\n';
    std::cout << l * u << '\n' << p * m2 << '\n';
    return 0;
}
