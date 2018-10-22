#pragma once

#include "coefficients.hpp"
#include "../utils/utils.hpp"
#include <functional>
#include <tuple>
#include <type_traits>

namespace ssh {

    namespace {

        template<typename Argument, typename Function, typename Value = Argument, typename... Parameters>
        struct caller {

            Function function;

            caller(const Function& function) : function(function) {}

            Value operator()(const Argument& argument, const Value& value, const std::tuple<Parameters...>& parameters) {
                return call(function, argument, value, parameters);
            }

            static Value call(const Function& function, const Argument& argument, const Value& value, const std::tuple<Parameters...>& parameters) {
                return std::apply([&function, &argument, &value](auto&&... parameters) { 
                    return function(argument, value, std::forward<decltype(parameters)>(parameters)...);
                }, parameters);
            }

        };
        
    }

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

    template<
        typename AVector, typename VVector, typename Function, typename Coeffs, 
        typename KVector = types::vector2d_t<double>,
        typename Argument = typename AVector::value_type, 
        typename Value = typename VVector::value_type,
        template<typename> typename PVector = types::vector1d_t, 
        typename KSubVector = typename KVector::value_type,
        typename... Parameters>
    class lazy_runge_kutta_p {

        public:

            lazy_runge_kutta_p() = delete;
            ~lazy_runge_kutta_p() = default;

            lazy_runge_kutta_p(const lazy_runge_kutta_p&) = default;
            lazy_runge_kutta_p(lazy_runge_kutta_p&&) = default;

            lazy_runge_kutta_p(const AVector& arguments, const VVector& initial_value, const Function& function,
                                const PVector<std::tuple<Parameters...>>& parameters, const Coeffs& coefficients) :
                _arguments(arguments), _parameters(parameters), _max_index(arguments.size() - 1), _initial_value(initial_value),
                _last_value(initial_value), _caller(_caller_t(function)), _coefficients(coefficients) {}

            lazy_runge_kutta_p(const AVector& arguments, const VVector& initial_value, const Function& function,
                                const PVector<std::tuple<Parameters...>>& parameters) :
                lazy_runge_kutta_p(arguments, initial_value, function, parameters, Coeffs()) {}

            VVector operator()(KVector& k) {
                if (!*this)
                    return _last_value;
                const auto arg = _arguments[_arg_index++];
                const auto d = _arguments[_arg_index] - arg;
                for (size_t i = 0; i < _parameters.size(); ++i) {
                    auto result = Value(0);
                    for (size_t j = 0; j < _coefficients.steps(); ++j) {
                        auto buff = Value(0);
                        for (size_t m = 0; m < j; ++m)
                            buff += _coefficients.get_a(j, m) * k[i][m];
                        result += _coefficients.get_b(j) * (
                            k[i][j] = _caller(arg + d * _coefficients.get_c(j), _last_value[i] + buff * d, _parameters[i]));
                    }
                    _last_value[i] += result * d;
                }
                return _last_value;
            }

            VVector operator()() {
                static auto k = std::is_constructible_v<KVector, size_t, KSubVector> 
                                    ? std::is_constructible_v<KSubVector, size_t> 
                                        ? KVector(_parameters.size(), KSubVector(_coefficients.steps()))
                                        : KVector(_parameters.size())
                                    : std::is_constructible_v<KVector, size_t>
                                        ? KVector(_parameters.size())
                                        : KVector();
                return (*this)(k);
            }

            operator bool() const {
                return _max_index > _arg_index;
            }

            void reset() {
                _arg_index = 0;
                _last_value = _initial_value;
            }

        private:

            typedef caller<Argument, Function, Value, Parameters...> _caller_t;

            const AVector _arguments;
            const PVector<std::tuple<Parameters...>> _parameters;
            size_t _arg_index = 0;
            const size_t _max_index;
            const AVector _initial_value;
            AVector _last_value;
            _caller_t _caller;
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
