#include <cmath>
#include <vector>
#include <fstream>
#include "dork.hpp"

//Solving initial value problem y' = f(x, y), y(a) = y0 using standard Runge Kutta method

int main() {
    //Computational area
    const double a = 0;
    const double b = 10;

    //Number of points
    const size_t n = 20;

    //Initial value
    const double y0 = 0;

    //Specify function
    const auto f = [](const double x, const double y) { return std::cos(x); };

    //Create uniform solver
    auto uniform_solver = dork::rk4<double>(a, b, n)(f)(y0);

    //Output solution
    auto out = std::ofstream("output_uniform.txt");
    for (const auto& [x, y] : uniform_solver)
        out << x << ' ' << y << '\n';

    //Use non uniformly distributed points
    const std::vector<double> points = { 0, 0.3, 1, 1.2, 2, 2.7, 3, 3.4, 4, 4.1, 5, 5.5, 6, 6.8, 7, 7.9, 8, 8.8, 9.6, 10 };

    //Create solver
    auto solver = dork::rk4<double>(points.begin(), points.end())(f)(y0);
    // or
    // auto solver = dork::rk4<double>(points)(f)(y0);

    //Output solution
    out = std::ofstream("output.txt");
    for (const auto& [x, y] : solver)
        out << x << ' ' << y << '\n';
}
