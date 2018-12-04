#include <cmath>
#include <vector>
#include <fstream>
#include "../include/runge_kutta/runge_kutta.hpp"
#include "../include/runge_kutta/coefficients.hpp"

//Solving initial value problem y' = f(x, y, params), y(a, params) = y0
//where params are some set of independent variables

int main() {
    //First create runge kutta method using some coefficients
    const auto method = dork::runge_kutta::new_method(dork::runge_kutta_coefficients<double, 4>());

    //Computational area
    const double a = 0;
    const double b = 10;

    //Number of points
    const size_t n = 20;

    //Create params for problem y' = f(x, y, u, v)
    const std::vector<std::tuple<double, double>> params = {
        std::make_tuple(M_PI, M_PI_2), // y' = f(x, y, M_PI, M_PI_2)
        std::make_tuple(M_PI_2, M_PI), // y' = f(x, y, M_PI_2, M_PI)
        std::make_tuple(0., 0.)        // y' = f(x, y, 0, 0)
    };

    //Specify function
    const auto f = [](const double x, const double y, const double u, const double v) { 
        return std::cos(x + u - v); 
    };

    //Calculate initial values
    const std::vector<double> y0 = {
        std::sin(a + std::get<0>(params[0]) - std::get<1>(params[0])),
        std::sin(a + std::get<0>(params[1]) - std::get<1>(params[1])),
        std::sin(a + std::get<0>(params[2]) - std::get<1>(params[2]))
    };

    //Solve the problem, distributing points uniformly in computational area
    const auto unifrom_solution = method.solve_uniform_p(a, b, n, y0, f, params);

    //Output the solution
    auto out = std::ofstream("output_uniform.txt");
    for (const auto& it : unifrom_solution) {
        for (const auto& y : it)
            out << y << ' ';
        out << '\n';
    }

    // Use non uniformly distributed points
    const std::vector<double> points = { 0, 0.3, 1, 1.2, 2, 2.7, 3, 3.4, 4, 4.1, 5, 5.5, 6, 6.8, 7, 7.9, 8, 8.8, 9.6, 10 };

    //Solve these points
    auto solution = method.solve_p(points, y0, f, params);

    //Output the solution
    out = std::ofstream("output.txt");
    for (const auto& it : solution) {
        for (const auto& y : it)
            out << y << ' ';
        out << '\n';
    }
}
