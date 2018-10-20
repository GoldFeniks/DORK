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
    const auto args = ssh::utils::create_mesh_bounds<ssh::types::vector1d_t<double>>(0., 10., 101);
    std::vector<std::array<double, 1>> params;
    for (size_t i = 0; i < 11; ++i)
        params.push_back({ i / 10. * M_PI_2 });
    std::vector<double> init;
    for (const auto& it : params)
        init.push_back(std::sin(it[0]));
    auto rk = ssh::lazy_runge_kutta_p(args, init, [](double x, double y, double z) { return std::cos(x + z); }, params, ssh::runge_kutta_coefficients<double, 4>());
    std::ofstream out("./source/cmake-build-debug/test.txt");
    for (const auto& it : init)
        out << it << ' ';
    out << '\n';
    while (rk) {
        for (const auto& it : rk()) {
            out << it << ' ';
        }
        out << '\n';
    }
    // auto rk = ssh::lazy_runge_kutta_uniform(0., 10., 21, 0., [](double x, double y) { return std::cos(x) + std::sin(x); }, ssh::runge_kutta_coefficients<double, 4>());
    // auto o = ssh::runge_kutta_dense_output::solve(rk, 0.5, 0., 11, ssh::runge_kutta_dense_output_coefficients<double, 3>());
    // std::ofstream out("test.txt");
    // for (const auto it : o)
        // out << it << ' ';
}



