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
    ssh::types::vector1d_t<std::tuple<double>> params;
    for (size_t i = 0; i < 11; ++i)
        params.push_back({ i / 10. * M_PI_2 });
    std::vector<double> init;
    for (const auto& it : params)
        init.push_back(std::sin(std::get<0>(it)));
    auto rk_method = ssh::runge_kutta(ssh::runge_kutta_coefficients<double, 4>());
    auto rk = rk_method.create_lazy_uniform_p(0., 10., 101, init, [](double x, double y, double z) { return std::cos(x + z); }, params);
    std::ofstream out("./test.txt");
    for (const auto& it : init)
        out << it << ' ';
    out << '\n';
    while (rk) {
        for (const auto& it : rk()) {
            out << it << ' ';
        }
        out << '\n';
    }
}



