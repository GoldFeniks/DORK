#include <cmath>
#include <vector>
#include <fstream>
#include "../include/runge_kutta/runge_kutta.hpp"
#include "../include/runge_kutta/coefficients.hpp"

//Solving initial value problem y' = f(x, y), y(a) = y0

int main() {
    //First create runge kutta method using some coefficients
    const auto method = dork::runge_kutta(dork::runge_kutta_coefficients<double, 4>());

    //Computational area
    const double a = 0;
    const double b = 10;

    //Number of points
    const size_t n = 20;

    //Initial value
    const double y0 = 0;

    //Specify function
    const auto f = [](const double x, const double y) { return std::cos(x); };

    //Solve the problem, distributing points uniformly in computational area
    const auto uniform_solution = method.solve_uniform(a, b, n, y0, f);

    //Output solution
    auto out = std::ofstream("output_uniform.txt");
    for (const auto& it : uniform_solution)
        out << it << '\n';

    //Use non uniformly distributed points
    const std::vector<double> points = { 0, 0.3, 1, 1.2, 2, 2.7, 3, 3.4, 4, 4.1, 5, 5.5, 6, 6.8, 7, 7.9, 8, 8.8, 9.6, 10 };

    //Solve at these points
    const auto solution = method.solve(points, y0, f);

    //Output solution
    out = std::ofstream("output.txt");
    for (const auto& it : solution)
        out << it << '\n';
}
