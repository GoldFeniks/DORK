# DORK
Dense Output for Runge Kutta

Project implements generalized Runge Kutta method with Dense Output support. The implementation is based on [Solving Ordinary Differential Equations I  by Ernst Hairer, Syvert P. Nørsett, Gerhard Wanner](https://www.springer.com/gp/book/9783540566700).

## Minimal example

```c++
#include <iostream>
#include "dork.hpp"

int main() {
    for (const auto& [x, y] : dork::rk4<double>(0, 1, 10)([](auto&...) { return 1; })(0.))
        std::cout << x << ' ' << y << '\n';
}
```
