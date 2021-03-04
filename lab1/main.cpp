#include <fstream>
#include <iostream>
#include "lib/matrix.h"

int main() {
    TMatrix m1, m2;
    std::ifstream is("../test.txt");
    is >> m1 >> m2;
    std::cout << m1 * m2 << '\n';
    return 0;
}
