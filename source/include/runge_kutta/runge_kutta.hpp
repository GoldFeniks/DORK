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

            Value operator()(const Argument& argument, const Value& value, const std::tuple<Parameters...>& parameters) const {
                return call(function, argument, value, parameters);
            }

            static Value call(const Function& function, const Argument& argument, const Value& value, const std::tuple<Parameters...>& parameters) {
                return std::apply([&function, &argument, &value](auto&&... parameters) { 
                    return function(argument, value, std::forward<decltype(parameters)>(parameters)...);
                }, parameters);
            }

        };

    }// namespace

    namespace implementation {

        template<typename Vector, typename Value, typename Function, typename Coeffs, typename KVector = types::vector1d_t<Value>>
        class lazy_runge_kutta {

        public:

            lazy_runge_kutta() = delete;
            ~lazy_runge_kutta() = default;

            lazy_runge_kutta(const lazy_runge_kutta&) = default;
            lazy_runge_kutta(lazy_runge_kutta&&) = default;

            lazy_runge_kutta(const Vector& arguments, const Value& initial_value,
                             const Function& function, const Coeffs& coefficients) :
                _arguments(arguments), _max_index(arguments.size() - 1), _initial_value(initial_value),
                _last_value(initial_value), _function(function), _coefficients(coefficients) {}

            lazy_runge_kutta(const Vector& arguments, const Value& initial_value, 
                             const Function& function) :
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
                static auto k = utils::constructor<KVector>::construct(
                    utils::arguments(_coefficients.steps()),
                    utils::no_arguments()
                );
                return (*this)(k);
            }

            operator bool() const {
                return _max_index > _arg_index;
            }

            void reset() {
                _last_value = _initial_value;
                _arg_index = 0;
            }

        private:

            const Vector _arguments;
            size_t _arg_index = 0;
            const size_t _max_index;
            const Value _initial_value;
            Value _last_value;
            const Function _function;
            const Coeffs _coefficients;

        };

        template<typename Argument, typename Function, typename Coeffs, typename Value = Argument, typename KVector = types::vector1d_t<Value>>
        class lazy_runge_kutta_uniform {

        public:

            lazy_runge_kutta_uniform() = delete;
            ~lazy_runge_kutta_uniform() = default;

            lazy_runge_kutta_uniform(const lazy_runge_kutta_uniform&) = default;
            lazy_runge_kutta_uniform(lazy_runge_kutta_uniform&&) = default;

            lazy_runge_kutta_uniform(const Argument& a, const Argument& b, const size_t n, const Value initial_value,
                                     const Function function, const Coeffs coefficients) :
                    _d((b - a) / (n - 1)), _left_argument(a), _right_argument(b), 
                    _last_argument(a), _initial_value(initial_value), _last_value(initial_value),
                    _function(function), _coefficients(coefficients) {}

            lazy_runge_kutta_uniform(const Argument& a, const Argument& b, const size_t n, const Value initial_value,
                                     const Function function) :
                    lazy_runge_kutta_uniform(a, b, n, initial_value, function, Coeffs()) {}

            Value operator()(KVector& k) {
                if (!*this)
                    return Value(0);
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
                static auto k = utils::constructor<KVector>::construct(
                    utils::arguments(_coefficients.steps()),
                    utils::no_arguments()
                );
                return (*this)(k);
            }

            operator bool() const {
                return _last_argument < _right_argument;
            }

        private:

            const Argument _d, _left_argument, _right_argument;
            Argument _last_argument;
            const Value _initial_value;
            Value _last_value;
            const Function _function;
            const Coeffs _coefficients;

        };

        template<
            typename AVector, typename VVector, typename Function, typename Coeffs, 
            typename Argument = typename AVector::value_type, 
            typename Value = typename VVector::value_type,
            typename KVector = types::vector2d_t<Value>,
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
                    static auto k = utils::constructor<KVector>::construct(
                        utils::arguments(_parameters.size(), 
                                        utils::constructor<KSubVector>::construct(
                                            utils::arguments(_coefficients.steps()), 
                                            utils::no_arguments())),
                        utils::arguments(_parameters.size()),
                        utils::no_arguments());
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
                const _caller_t _caller;
                const Coeffs _coefficients;

        }; 

        template<
            typename Argument, typename Function, typename Coeffs, 
            typename VVector = types::vector1d_t<Argument>, 
            typename Value = typename VVector::value_type,
            typename KVector = types::vector2d_t<Value>,
            template<typename> typename PVector = types::vector1d_t, 
            typename KSubVector = typename KVector::value_type,
            typename... Parameters>
        class lazy_runge_kutta_uniform_p {

            public:

                lazy_runge_kutta_uniform_p() = delete;
                ~lazy_runge_kutta_uniform_p() = default;

                lazy_runge_kutta_uniform_p(const lazy_runge_kutta_uniform_p&) = default;
                lazy_runge_kutta_uniform_p(lazy_runge_kutta_uniform_p&&) = default;

                lazy_runge_kutta_uniform_p(const Argument& a, const Argument& b, const size_t n, 
                                           const VVector& initial_value, const Function& function, 
                                           const PVector<std::tuple<Parameters...>>& parameters, 
                                           const Coeffs& coefficients) :
                    _d((b - a) / Argument(n)), _left_argument(a), _right_argument(b), _last_argument(a), 
                    _initial_value(initial_value), _last_value(initial_value), _parameters(parameters),
                    _caller(_caller_t(function)), _coefficients(coefficients) {}

                lazy_runge_kutta_uniform_p(const Argument& a, const Argument& b, const size_t n, 
                                           const VVector& initial_value, const Function& function, 
                                           const PVector<std::tuple<Parameters...>>& parameters) :
                    lazy_runge_kutta_uniform_p(a, b, n, initial_value, function, parameters, Coeffs()) {}

                VVector operator()(KVector& k) {
                    if (!*this)
                        return _last_value;
                    for (size_t i = 0; i < _parameters.size(); ++i) {
                        auto result = Value(0);
                        for (size_t j = 0; j < _coefficients.steps(); ++j) {
                            auto buff = Value(0);
                            for (size_t m = 0; m < j; ++m)
                                buff += _coefficients.get_a(j, m) * k[i][m];
                            result += _coefficients.get_b(j) * (
                                k[i][j] = _caller(_last_argument + _d * _coefficients.get_c(j), _last_value[i] + buff * _d, _parameters[i]));
                        }
                        _last_value[i] += result * _d;
                    }
                    _last_argument += _d;
                    return _last_value;
                }

                VVector operator()() {
                    static auto k = utils::constructor<KVector>::construct(
                        utils::arguments(_parameters.size(), 
                                        utils::constructor<KSubVector>::construct(
                                            utils::arguments(_coefficients.steps()), 
                                            utils::no_arguments())),
                        utils::arguments(_parameters.size()),
                        utils::no_arguments());
                    return (*this)(k);
                }

                operator bool() const {
                    return _last_argument < _right_argument;
                }

                void reset() {
                    _last_argument = _left_argument;
                    _last_value = _initial_value;
                }

            private:

                typedef caller<Argument, Function, Value, Parameters...> _caller_t;

                const Argument _d, _left_argument, _right_argument;
                Argument _last_argument;
                const VVector _initial_value;
                VVector _last_value;
                const PVector<std::tuple<Parameters...>> _parameters;
                const _caller_t _caller;
                const Coeffs _coefficients;

        };
        
    }// namespace implementation

    template<typename Coeffs, typename Argument = typename Coeffs::value_type, typename Value = Argument>
    class runge_kutta {

        public:

            runge_kutta() : _coefficients(Coeffs()) {}
            runge_kutta(const Coeffs& coefficients) : _coefficients(coefficients) {}

            template<typename Function, typename Vector = types::vector1d_t<Argument>, typename KVector = types::vector1d_t<Value>>
            auto create_lazy(const Vector& arguments, const Value& initial_value, const Function& function) {
                return implementation::
                    lazy_runge_kutta<Vector, Value, Function, Coeffs, KVector>(arguments, initial_value, function, _coefficients);
            }

            template<typename Function, typename KVector = types::vector1d_t<Value>>
            auto create_lazy_uniform(const Argument& a, const Argument& b, size_t n, const Value& initial_value, const Function& function) {
                return implementation::
                    lazy_runge_kutta_uniform<Argument, Function, Coeffs, Value, KVector>(a, b, n, initial_value, function, _coefficients);
            }

            template<
                typename Function,
                typename KVector = types::vector2d_t<Value>,
                typename AVector = types::vector1d_t<Argument>,
                typename VVector = types::vector1d_t<Value>,
                template<typename> typename PVector = types::vector1d_t, 
                typename KSubVector = typename KVector::value_type,
                typename... Parameters>
            auto create_lazy_p(const AVector& arguments, const VVector& initial_value, const Function& function,
                               const PVector<std::tuple<Parameters...>>& parameters) 
            {
                return implementation::
                    lazy_runge_kutta_p<AVector, VVector, Function, Coeffs, Argument, Value, KVector, PVector, KSubVector, Parameters...>
                        (arguments, initial_value, function, parameters, _coefficients);
            }

            template<
                typename Function,
                typename KVector = types::vector2d_t<Value>,
                typename VVector = types::vector1d_t<Argument>, 
                template<typename> typename PVector = types::vector1d_t, 
                typename KSubVector = typename KVector::value_type,
                typename... Parameters>
            auto create_lazy_uniform_p(const Argument& a, const Argument& b, const size_t n, const VVector& initial_value,
                                       const Function& function, const PVector<std::tuple<Parameters...>>& parameters) 
            {
                return implementation::
                    lazy_runge_kutta_uniform_p<Argument, Function, Coeffs, VVector, Value, KVector, PVector, KSubVector, Parameters...>
                        (a, b, n, initial_value, function, parameters, _coefficients); 
            }

        private:

            const Coeffs _coefficients;

    };

}// namespace ssh
