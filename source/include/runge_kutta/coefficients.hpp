#pragma once

#include "../utils/types.hpp"
#include <cstddef>

namespace ssh {

    template<typename Vector2, typename Vector, typename Value = typename Vector::value_type, typename Index = size_t>
    class coefficients {

    public:

        coefficients() = default;
        ~coefficients() = default;

        coefficients(const Vector2& a_table, const Vector& b_table, const Vector& c_table)
            : _a_table(a_table), _b_table(b_table), _c_table(c_table) {}

        const Value& get_a(const Index i, const Index j) const {
            return _a_table[i][j];
        }

        const Value& get_b(const Index i) const {
            return _b_table[i];
        }

        const Value get_c(const Index i) const {
            return _c_table[i];
        }

        const Index steps() const {
            return _a_table.size();
        }

    private:

        const Vector2 _a_table;
        const Vector _b_table, _c_table;

    };

    template<typename T, size_t N>
    class runge_kutta_coefficients : public coefficients<types::array2d_t<T, N>, types::array1d_t<T, N>> {};

    template<typename T>
    class runge_kutta_coefficients<T, 1> : public coefficients<types::array2d_t<T, 1>, types::array1d_t<T, 1>> {

    public:

        static constexpr types::array2d_t<T, 1> a_table = { T(1) };
        static constexpr types::array1d_t<T, 1> b_table = { T(1) };
        static constexpr types::array1d_t<T, 1> c_table = { T(0) };

        runge_kutta_coefficients() : coefficients<types::array2d_t<T, 1>, types::array1d_t<T, 1>>(a_table, b_table, c_table) {};

    };

    template<typename T>
    class runge_kutta_coefficients<T, 2> : public coefficients<types::array2d_t<T, 2>, types::array1d_t<T, 2>> {

    public:

        static constexpr types::array2d_t<T, 2> a_table = { T(0), T(0),  T(1) / T(2), T(0)  };
        static constexpr types::array1d_t<T, 2> b_table = { T(0), T(1) };
        static constexpr types::array1d_t<T, 2> c_table = { T(0), T(1) / T(2) };

        runge_kutta_coefficients() : coefficients<types::array2d_t<T, 2>, types::array1d_t<T, 2>>(a_table, b_table, c_table) {};

    };

    template<typename T>
    class runge_kutta_coefficients<T, 4> : public coefficients<types::array2d_t<T, 4>, types::array1d_t<T, 4>> {

    public:

        static constexpr types::array2d_t<T, 4> a_table = {
                T(0),        T(0),        T(0), T(0),
                T(1) / T(2), T(0),        T(0), T(0),
                T(0),        T(1) / T(2), T(0), T(0),
                T(0),        T(0),        T(1), T(0)
        };
        static constexpr types::array1d_t<T, 4> b_table = {
                T(1) / T(6),
                T(1) / T(3),
                T(1) / T(3),
                T(1) / T(6)
        };
        static constexpr types::array1d_t<T, 4> c_table = {
                T(0),
                T(1) / T(2),
                T(1) / T(2),
                T(1)
        };

        runge_kutta_coefficients() : coefficients<types::array2d_t<T, 4>, types::array1d_t<T, 4>>(a_table, b_table, c_table) {}
    };


}// namespace ssh
