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
    const auto rk_method = ssh::runge_kutta(ssh::runge_kutta_coefficients<double, 4>());
    const auto res1 = rk_method.solve_uniform(0., 10., 100,  0., [](const double& x, const double& y) { return  std::cos(x); });
    const auto res2 = rk_method.solve_uniform(0., 10., 100,  1., [](const double& x, const double& y) { return -std::sin(x); });
    const auto res3 = rk_method.solve_uniform(0., 10., 100,  0., [](const double& x, const double& y) { return -std::cos(x); });
    const auto res4 = rk_method.solve_uniform(0., 10., 100, -1., [](const double& x, const double& y) { return  std::sin(x); });
    std::ofstream out("./source/cmake-build-debug/test.txt");
    for (size_t i = 0; i < 100; ++i)
        out << res1[i] << ' ' << res2[i] << ' ' << res3[i] << ' ' << res4[i] << '\n';
}
