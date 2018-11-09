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

    //Create lazy solver, distributing points uniformly in computational area
    auto solver_uniform = method.create_lazy_uniform(a, b, n, y0, f);

    //Calculate and output the solution
    auto out = std::ofstream("output_uniform.txt");
    out << y0 << '\n'; //Solver does not output initial value
    while (solver_uniform)
        out << solver_uniform() << '\n';

    //Use non uniformly distributed points
    const std::vector<double> points = { 0, 0.3, 1, 1.2, 2, 2.7, 3, 3.4, 4, 4.1, 5, 5.5, 6, 6.8, 7, 7.9, 8, 8.8, 9.6, 10 };

    //Create solver for these points
    auto solver = method.create_lazy(points, y0, f);

    //Calculate and output the solution
    out = std::ofstream("output.txt");
    out << y0 << '\n'; //Solver does not output initial value
    while (solver)
        out << solver() << '\n';
}
