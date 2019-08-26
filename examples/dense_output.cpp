#include <cmath>
#include <vector>
#include <fstream>
#include "dork.hpp"

//Solving initial value problem y' = f(x, y), y(a) = y0 using DOPRI method and dense output

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

    //Create solver with dense output support. If method does not support dense output linear interpolation is used instead
    auto solver = dork::dopri5_do5<double>(a, b, n)(f)(y0);

    //Output solution
    auto out = std::ofstream("output.txt");
    while (solver) {
        solver.next();
        //Output 10 values in each step
        for (size_t i = 0; i < 10; ++i) {
            const auto [x, y] = solver(0.1 * i);
            out << x << ' ' << y << '\n';
        }
    }
    //Output last value
    const auto [x, y] = solver(1);
    out << x << ' ' << y << '\n';
}
