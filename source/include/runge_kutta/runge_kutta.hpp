#pragma once

#include "coefficients.hpp"
#include "../utils/utils.hpp"
#include <functional>

namespace ssh {

    template<typename Value, typename Vector, typename Function, typename Coeffs, typename KVector = types::vector1d_t<Value>>
    class lazy_runge_kutta {

    public:

        lazy_runge_kutta(const Vector arguments, const Value initial_value,
                         const Function function, const Coeffs coefficients) :
            _arguments(arguments), _max_index(arguments.size() - 1), _last_value(initial_value),
            _function(function), _coefficients(coefficients) {}

        lazy_runge_kutta(const Vector arguments, const Value initial_value, const Function function) :
            lazy_runge_kutta(arguments, initial_value, function, Coeffs()) {}

        Value operator()(KVector& k) {
            if (!*this)
                return Value(0);
            auto result = Value(0);
            const auto arg = _arguments[_arg_index++];
            const auto d = _arguments[_arg_index] - arg;
            for (size_t i = 0; i < _coefficients.steps(); ++i) {
                auto buff = Value(0);
                for (size_t j = 0; j < i; ++j)
                    buff += _coefficients.get_a(i, j) * k[j];
                result += _coefficients.get_b(i) * (k[i] = _function(arg + d * _coefficients.get_c(i), _last_value + buff * d));
            }
            return _last_value += result * d;
        }

        Value operator()() {
            static KVector k(_coefficients.steps());
            return (*this)(k);
        }

        operator bool() const {
            return _max_index > _arg_index;
        }

    private:

        const Vector _arguments;
        size_t _arg_index = 0;
        const size_t _max_index;
        Value _last_value;
        const Function _function;
        const Coeffs _coefficients;

    };

    template<typename Value, typename Argument, typename Function, typename Coeffs, typename KVector = types::vector1d_t<Value>>
    class lazy_runge_kutta_uniform {

    public:

        lazy_runge_kutta_uniform(const Argument& a, const Argument& b, const size_t n, const Value initial_value,
                                 const Function function, const Coeffs coefficients) :
                _d((b - a) / (n - 1)), _last_argument(a), _max_index(n), _last_value(initial_value),
                _function(function), _coefficients(coefficients) {}

        lazy_runge_kutta_uniform(const Argument& a, const Argument& b, const size_t n, const Value initial_value,
                         const Function function) :
                lazy_runge_kutta_uniform(a, b, n, initial_value, function, Coeffs()) {}

        Value operator()(KVector& k) {
            if (!*this)
                return Value(0);
            ++_index;
            auto result = Value(0);
            for (size_t i = 0; i < _coefficients.steps(); ++i) {
                auto buff = Value(0);
                for (size_t j = 0; j < i; ++j)
                    buff += _coefficients.get_a(i, j) * k[j];
                result += _coefficients.get_b(i) *
                        (k[i] = _function(_last_argument + _d * _coefficients.get_c(i), _last_value + buff * _d));
            }
            _last_argument += _d;
            return _last_value += result * _d;
        }

        Value operator()() {
            static KVector k(_coefficients.steps());
            return (*this)(k);
        }

        operator bool() const {
            return _max_index > _index;
        }

    private:

        const Argument _d;
        Argument _last_argument;
        size_t _index = 1;
        const size_t _max_index;
        Value _last_value;
        const Function _function;
        const Coeffs _coefficients;

    };

    template<typename Value, typename Argument, typename Coeffs>
    class runge_kutta {

    public:

        template<typename Function, typename Vector = types::vector1d_t<Value>>
        static Vector solve(const Argument& a, const Argument& b,  const size_t n,
                            const Value& initial_value, const Function& function, const Coeffs& coefficients) {
            Vector result(n);
            solve(result, a, b, n, initial_value, function, coefficients);
            return result;
        }

        template<typename Function, typename Vector = types::vector1d_t<Value>>
        static void solve(Vector& vector, const Argument& a, const Argument& b, const size_t n,
                          const Value& initial_value, const Function& function, const Coeffs& coefficients) {
            vector[0] = initial_value;
            auto rk = lazy_runge_kutta_uniform(a, b, n, initial_value, function, coefficients);
            for (size_t i = 1; i < n; ++i)
                vector[i] = rk();
        }

        template<typename Vector, typename Function, typename RVector = types::vector1d_t<Value>>
        static RVector solve(const Vector& arguments, const Value& initial_value,
                             const Function& function, const Coeffs& coefficients) {
            RVector result(arguments.size());
            solve(result, arguments, initial_value, function, coefficients);
            return result;
        }

        template<typename Vector, typename Function, typename RVector = types::vector1d_t<Value>>
        static void solve(RVector& vector, const Vector& arguments, const Value& initial_value,
                             const Function& function, const Coeffs& coefficients) {
            vector[0] = initial_value;
            auto rk = lazy_runge_kutta(arguments, initial_value, function, coefficients);
            for (size_t i = 1; i < arguments.size(); ++i)
                vector[i] = rk();
        }

    };

}// namespace ssh
