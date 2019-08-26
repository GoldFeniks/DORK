#include <cmath>
#include <vector>
#include <fstream>
#include "dork.hpp"

//Solving initial value problem system
//    y' = f(x, y, z), y(a) = y0,
//    z' = g(x, y, z), z(a) = z0
// using standard Runke Kutta method

int main() {
    //Computational area
    const double a = 0;
    const double b = 10;

    //Number of points
    const size_t n = 20;

    //Initial values
    const double y0 =  1;
    const double z0 = -1;

    //Specify functions
    const auto f = [](const double& x, const double& y, const double& z) { return std::sin(x) * z; };
    const auto g = [](const double& x, const double& y, const double& z) { return std::sin(x) * y; };

    //Create solver
    auto solver = dork::rk4<double>(a, b, n)(f, g)(y0, z0);

    //Output solution
    auto out = std::ofstream("output.txt");
    for (const auto& [x, y, z] : solver)
        out << x << ' ' << y << ' ' << z << '\n';
}
