#include <cmath>
#include <vector>
#include <fstream>
#include "dork.hpp"

//Solving initial value problem y' = f(x, y), y(a) = y0 using DOPRI method with step control

int main() {
    //Computational area
    const double a = 0;
    const double b = 10;

    //Inital step
    const double h = 0.1;

    //Specify tolerances. Smaller values correspond to smaller step size
    double atol = 1e-5; //Absolute error tolerance
    double rtol = 1e-5; //Reference error tolerance

    //Initial value
    const double y0 = 0;

    //Specify function
    const auto f = [](const double x, const double y) { return std::cos(x); };

    //Create solver
    auto solver1 = dork::dopri5<double>(a, b, h, atol, rtol)(f)(y0);

    //Output solution
    auto out = std::ofstream("output1.txt");
    for (const auto& [x, y] : solver1)
        out << x << ' ' << y << '\n';

    //Solve with smaller tolerances
    atol = 1e-10;
    rtol = 1e-10;
    auto solver2 = dork::dopri5<double>(a, b, h, atol, rtol)(f)(y0);

    out = std::ofstream("output2.txt");
    for (const auto& [x, y] : solver2)
        out << x << ' ' << y << '\n';
}
