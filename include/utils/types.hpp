#pragma once

#include <vector>
#include <array>

namespace dork {

    namespace types {

        namespace implementation {

            template<typename T, size_t N>
            struct vector_type {

                using type = std::vector<typename vector_type<T, N - 1>::type>;

            };

            template<typename T>
            struct vector_type<T, 0> {

                using type = T;

            };

            template<typename, size_t...>
            struct array_type;

            template<typename T, size_t N, size_t... Tail>
            struct array_type<T, N, Tail...> {

                using type = std::array<typename array_type<T, Tail...>::type, N>;

            };

            template<typename T>
            struct array_type<T> {

                using type = T;

            };

        };

        template<typename T, size_t N>
        using vector_t = typename implementation::vector_type<T, N>::type;

        template<typename T, size_t... N>
        using array_t = typename implementation::array_type<T, N...>::type;

        template<typename T>
        using vector1d_t = std::vector<T>;

        template<typename T>
        using vector2d_t = vector_t<T, 2>;

        template<typename T, size_t N>
        using array1d_t = std::array<T, N>;

        template<typename T, size_t N>
        using array2d_t = array_t<T, N, N>;

    }

}// namespace dork
