#pragma once

#include "../utils/types.hpp"
#include <cstddef>
#include <numeric>

namespace ssh {

    template<typename Vector2, typename Vector, typename Value = typename Vector::value_type, typename Index = size_t>
    class coefficients {

    public:

        using value_type = Value;

        coefficients() = default;
        ~coefficients() = default;

        coefficients(const Vector2& a_table, const Vector& b_table, const Vector& c_table)
            : _a_table(a_table), _b_table(b_table), _c_table(c_table) {}

        coefficients(const Vector2& a_table, const Vector& b_table)
            : _a_table(a_table), _b_table(b_table), _c_table(_get_c_table(a_table, _b_table)) {}

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

        static Vector _get_c_table(const Vector2& a_table, const Vector& b_table) {
            Vector result = b_table;
            for (Index i = 0; i < result.size(); ++i)
                result[i] = std::accumulate(a_table[i].begin(), a_table[i].end(), Value(0));
            return result;
        }

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

        runge_kutta_coefficients() : coefficients<types::array2d_t<T, 2>, types::array1d_t<T, 2>>(a_table, b_table) {};

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

        runge_kutta_coefficients() : coefficients<types::array2d_t<T, 4>, types::array1d_t<T, 4>>(a_table, b_table) {}
    };

    template<typename T>
    class dormand_prince_coefficients1 : public coefficients<types::array2d_t<T, 7>, types::array1d_t<T, 7>> {

    public:

        static constexpr types::array2d_t<T, 7> table =
                {
                        T(0),                T(0),               T(0),                T(0),             T(0),               T(0),          T(0),
                        T(1)     / T(5),     T(0),               T(0),                T(0),             T(0),               T(0),          T(0),
                        T(3)     / T(40),    T(9)     / T(40),   T(0),                T(0),             T(0),               T(0),          T(0),
                        T(44)    / T(45),   -T(56)    / T(15),   T(32)    / T(9),     T(0),             T(0),               T(0),          T(0),
                        T(19372) / T(6561), -T(25360) / T(2187), T(64448) / T(6561), -T(212) / T(729),  T(0),               T(0),          T(0),
                        T(9017)  / T(3168), -T(355)   / T(33),   T(46732) / T(5247),  T(49)  / T(176), -T(5103) / T(18656), T(0),          T(0),
                        T(35)    / T(384),   T(0),               T(500)   / T(1113),  T(125) / T(192), -T(2187) / T(6784),  T(11) / T(84), T(0),
                };

        dormand_prince_coefficients1() : coefficients<types::array2d_t<T, 7>, types::array1d_t<T, 7>>(table, table[6]) {}
    };

    template<typename T>
    class dormand_prince_coefficients2 : public coefficients<types::array2d_t<T, 7>, types::array1d_t<T, 7>> {

    public:

        static constexpr types::array2d_t<T, 7> a_table = dormand_prince_coefficients1<T>::table;
        static constexpr types::array1d_t<T, 7> b_table =
                {
                         T(5179)  / T(57600),
                         T(0),
                         T(7571)  / T(16695),
                         T(393)   / T(640),
                        -T(92097) / T(339200),
                         T(187)   / T(2100),
                         T(1)     / T(40)
                };

        dormand_prince_coefficients2() : coefficients<types::array2d_t<T, 7>, types::array1d_t<T, 7>>(a_table, b_table) {}

    };

}// namespace ssh
