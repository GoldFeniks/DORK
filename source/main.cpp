#include "include/utils/utils.hpp"
#include "include/runge_kutta/coefficients.hpp"
#include "include/runge_kutta/runge_kutta.hpp"
#include "include/runge_kutta/dense_output.hpp"
#include "include/runge_kutta/dense_output_coefficients.hpp"
#include "include/utils/constructor.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

int main() {
    const auto rk_method = ssh::runge_kutta(ssh::dormand_prince_coefficients1<double>());
    const auto deop = ssh::dense_output::from_coefficients(ssh::dormand_prince_dense_output_coefficients1<double>());
    ssh::types::vector1d_t<std::tuple<double>> params;
    for (size_t i = 0; i < 11; ++i)
        params.push_back({ i / 10. * M_PI_2 });
    std::vector<double> init;
    for (const auto& it : params)
        init.push_back(std::sin(std::get<0>(it)));
    auto rk = rk_method.create_lazy_uniform_p(0., 10., 11, init, [](double x, double y, double z) { return std::cos(x + z); }, params);
    std::ofstream out("./source/cmake-build-debug/test.txt");
    const auto res = deop.solve_uniform_p(rk, 1000);
    for (const auto& it : res) {
        for (const auto& itt : it) 
            out << itt << ' ';
        out << '\n';
    }
}
