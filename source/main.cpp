#include "include/utils/utils.hpp"
#include "include/runge_kutta/coefficients.hpp"
#include "include/runge_kutta/runge_kutta.hpp"
#include "include/runge_kutta/dense_output.hpp"
#include "include/runge_kutta/dense_output_coefficients.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

int main() {
    auto rk = ssh::lazy_runge_kutta_uniform(0., 10., 21, 0., [](double x, double y) { return std::cos(x) + std::sin(x); }, ssh::runge_kutta_coefficients<double, 4>());
    auto o = ssh::runge_kutta_dense_output::solve(rk, 0.5, 0., 11, ssh::runge_kutta_dense_output_coefficients<double, 3>());
    std::ofstream out("test.txt");
    for (const auto it : o)
        out << it << ' ';
}



