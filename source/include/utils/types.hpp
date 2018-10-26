#pragma once

#include <vector>
#include <array>

namespace ssh {

    namespace types {

        template<typename T>
        using vector1d_t = std::vector<T>;

        template<typename T>
        using vector2d_t = std::vector<std::vector<T>>;

        template<typename T, size_t N>
        using array1d_t = std::array<T, N>;

        template<typename T, size_t N>
        using array2d_t = std::array<std::array<T, N>, N>;

    }

}// namespace ssh
