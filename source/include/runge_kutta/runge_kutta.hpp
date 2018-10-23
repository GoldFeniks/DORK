#pragma once

#include "coefficients.hpp"
#include "../utils/utils.hpp"
#include <tuple>
#include <type_traits>

namespace ssh {

    namespace {

        template<typename Argument, typename Function, bool  StoreRef = false, typename Value = Argument, typename... Parameters>
        struct caller {

            std::conditional_t<StoreRef, const Function&, const Function> function;

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

        template<
            typename Vector, typename Value, typename Function, typename Coeffs,
            bool StoreRef = false, typename KVector = types::vector1d_t<Value>>
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
                lazy_runge_kutta(arguments, initial_value, function, Coeffs()) 
            {
                static_assert(!StoreRef, "Coeffitients must be provided to store reference");
            }

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

            auto size() const {
                return _arguments.size();
            }

            auto initial_value() const {
                return _initial_value;
            }

        private:

            std::conditional_t<StoreRef, const Vector&, const Vector> _arguments;
            size_t _arg_index = 0;
            const size_t _max_index;
            std::conditional_t<StoreRef, const Value&, const Value> _initial_value;
            Value _last_value;
            std::conditional_t<StoreRef, const Function&, const Function> _function;
            std::conditional_t<StoreRef, const Coeffs&, const Coeffs> _coefficients;

        };

        template<
            typename Argument, typename Function, typename Coeffs, 
            bool StoreRef = false,
            typename Value = Argument, 
            typename KVector = types::vector1d_t<Value>>
        class lazy_runge_kutta_uniform {

        public:

            lazy_runge_kutta_uniform() = delete;
            ~lazy_runge_kutta_uniform() = default;

            lazy_runge_kutta_uniform(const lazy_runge_kutta_uniform&) = default;
            lazy_runge_kutta_uniform(lazy_runge_kutta_uniform&&) = default;

            lazy_runge_kutta_uniform(const Argument& a, const Argument& b, const size_t n, const Value& initial_value,
                                     const Function& function, const Coeffs& coefficients) :
                    _d((b - a) / (n - 1)), _n(n), _left_argument(a), _right_argument(b), 
                    _last_argument(a), _initial_value(initial_value), _last_value(initial_value),
                    _function(function), _coefficients(coefficients) {}

            lazy_runge_kutta_uniform(const Argument& a, const Argument& b, const size_t n, const Value& initial_value,
                                     const Function& function) :
                    lazy_runge_kutta_uniform(a, b, n, initial_value, function, Coeffs()) 
            {
                static_assert(!StoreRef, "Coeffitients must be provided to store reference");
            }

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

            void reset() {
                _last_argument = _left_argument;
                _last_value = _initial_value;
            }

            auto size() const {
                return _n;
            }

            auto initial_value() const {
                return _initial_value;
            }

        private:

            const Argument _d;
            const size_t _n;
            std::conditional_t<StoreRef, const Argument&, const Argument> _left_argument, _right_argument;
            Argument _last_argument;
            std::conditional_t<StoreRef, const Value&, const Value> _initial_value;
            Value _last_value;
            std::conditional_t<StoreRef, const Function&, const Function> _function;
            std::conditional_t<StoreRef, const Coeffs&, const Coeffs> _coefficients;

        };

        template<
            typename AVector, typename VVector, typename Function, typename Coeffs, 
            bool StoreRef = false,
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
                    lazy_runge_kutta_p(arguments, initial_value, function, parameters, Coeffs()) 
                {
                    static_assert(!StoreRef, "Coeffitients must be provided to store reference");
                }

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

                auto size() const {
                    return _arguments.size();
                }

                auto initial_value() const {
                    return _initial_value;
                }

            private:

                typedef caller<Argument, Function, StoreRef, Value, Parameters...> _caller_t;

                std::conditional_t<StoreRef, const AVector&, const AVector> _arguments;
                std::conditional_t<StoreRef, 
                    const PVector<std::tuple<Parameters...>>&, 
                    const PVector<std::tuple<Parameters...>>> _parameters;
                size_t _arg_index = 0;
                const size_t _max_index;
                std::conditional_t<StoreRef, const AVector&, const AVector> _initial_value;
                AVector _last_value;
                const _caller_t _caller;
                std::conditional_t<StoreRef, const Coeffs&, const Coeffs> _coefficients;

        }; 

        template<
            typename Argument, typename Function, typename Coeffs, 
            bool StoreRef = false,
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
                    _d((b - a) / Argument(n - 1)), _n(n), _left_argument(a), _right_argument(b), _last_argument(a), 
                    _initial_value(initial_value), _last_value(initial_value), _parameters(parameters),
                    _caller(_caller_t(function)), _coefficients(coefficients) {}

                lazy_runge_kutta_uniform_p(const Argument& a, const Argument& b, const size_t n, 
                                           const VVector& initial_value, const Function& function, 
                                           const PVector<std::tuple<Parameters...>>& parameters) :
                    lazy_runge_kutta_uniform_p(a, b, n, initial_value, function, parameters, Coeffs()) 
                {
                    static_assert(!StoreRef, "Coeffitients must be provided to store reference");
                }

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

                auto size() const {
                    return _n;
                }

                auto initial_value() const {
                    return _initial_value;
                }

            private:

                typedef caller<Argument, Function, StoreRef, Value, Parameters...> _caller_t;

                const Argument _d;
                const size_t _n;
                std::conditional_t<StoreRef, const Argument&, const Argument> _left_argument, _right_argument;
                Argument _last_argument;
                std::conditional_t<StoreRef, const VVector&, const VVector> _initial_value;
                VVector _last_value;
                std::conditional_t<StoreRef, 
                    const PVector<std::tuple<Parameters...>>&,
                    const PVector<std::tuple<Parameters...>>> _parameters;
                const _caller_t _caller;
                std::conditional_t<StoreRef, const Coeffs&, const Coeffs> _coefficients;

        };
        
    }// namespace implementation

    template<typename Coeffs, typename Argument = typename Coeffs::value_type, typename Value = Argument>
    class runge_kutta {

        public:

            runge_kutta() : _coefficients(Coeffs()) {}
            runge_kutta(const Coeffs& coefficients) : _coefficients(coefficients) {}

            template<typename Function, typename Vector = types::vector1d_t<Argument>, typename KVector = types::vector1d_t<Value>>
            auto create_lazy(const Vector& arguments, const Value& initial_value, const Function& function) const {
                return implementation::
                    lazy_runge_kutta<Vector, Value, Function, Coeffs, false, KVector>(arguments, initial_value, function, _coefficients);
            }

            template<typename Function, typename KVector = types::vector1d_t<Value>>
            auto create_lazy_uniform(const Argument& a, const Argument& b, size_t n, const Value& initial_value, const Function& function) const {
                return implementation::
                    lazy_runge_kutta_uniform<Argument, Function, Coeffs, false, Value, KVector>(a, b, n, initial_value, function, _coefficients);
            }

            template<typename Function, typename Vector = types::vector1d_t<Argument>, typename KVector = types::vector1d_t<Value>>
            auto create_lazy_bound(const Vector& arguments, const Value& initial_value, const Function& function) const {
                return implementation::
                    lazy_runge_kutta<Vector, Value, Function, Coeffs, true, KVector>(arguments, initial_value, function, _coefficients);
            }

            template<typename Function, typename KVector = types::vector1d_t<Value>>
            auto create_lazy_uniform_bound(const Argument& a, const Argument& b, size_t n, const Value& initial_value, const Function& function) const {
                return implementation::
                    lazy_runge_kutta_uniform<Argument, Function, Coeffs, true, Value, KVector>(a, b, n, initial_value, function, _coefficients);
            }

            template<
                typename Function, 
                typename Vector = types::vector1d_t<Argument>, 
                typename RVector = types::vector1d_t<Value>, 
                typename KVector = types::vector1d_t<Value>>
            auto solve(const Vector& arguments, const Value& initial_value, const Function& function) const {
                return solver<RVector>::solve(create_lazy_bound<Function, Vector, KVector>(arguments, initial_value, function));
            }

            template<
                typename Function, 
                typename RVector = types::vector1d_t<Value>, 
                typename KVector = types::vector1d_t<Value>>
            auto solve_uniform(const Argument& a, const Argument& b, const size_t n, const Value& initial_value, const Function& function) const {
                return solver<RVector>::solve(create_lazy_uniform_bound<Function, KVector>(a, b, n, initial_value, function));
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
                               const PVector<std::tuple<Parameters...>>& parameters) const
            {
                return implementation::
                    lazy_runge_kutta_p<AVector, VVector, Function, Coeffs, false, Argument, Value, KVector, PVector, KSubVector, Parameters...>
                        (arguments, initial_value, function, parameters, _coefficients);
            }

            template<
                typename Function,
                typename KVector = types::vector2d_t<Value>,
                typename AVector = types::vector1d_t<Argument>,
                typename VVector = types::vector1d_t<Value>,
                template<typename> typename PVector = types::vector1d_t, 
                typename KSubVector = typename KVector::value_type,
                typename... Parameters>
            auto create_lazy_bound_p(const AVector& arguments, const VVector& initial_value, const Function& function,
                               const PVector<std::tuple<Parameters...>>& parameters) const
            {
                return implementation::
                    lazy_runge_kutta_p<AVector, VVector, Function, Coeffs, true, Argument, Value, KVector, PVector, KSubVector, Parameters...>
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
                                       const Function& function, const PVector<std::tuple<Parameters...>>& parameters) const
            {
                return implementation::
                    lazy_runge_kutta_uniform_p<Argument, Function, Coeffs, false, VVector, Value, KVector, PVector, KSubVector, Parameters...>
                        (a, b, n, initial_value, function, parameters, _coefficients); 
            }

            template<
                typename Function,
                typename KVector = types::vector2d_t<Value>,
                typename VVector = types::vector1d_t<Argument>, 
                template<typename> typename PVector = types::vector1d_t, 
                typename KSubVector = typename KVector::value_type,
                typename... Parameters>
            auto create_lazy_uniform_bound_p(const Argument& a, const Argument& b, const size_t n, const VVector& initial_value,
                                             const Function& function, const PVector<std::tuple<Parameters...>>& parameters) const
            {
                return implementation::
                    lazy_runge_kutta_uniform_p<Argument, Function, Coeffs, true, VVector, Value, KVector, PVector, KSubVector, Parameters...>
                        (a, b, n, initial_value, function, parameters, _coefficients); 
            }

            template<
                typename Function,
                typename RVector = types::vector2d_t<Value>,
                typename KVector = types::vector2d_t<Value>,
                typename AVector = types::vector1d_t<Argument>,
                typename VVector = types::vector1d_t<Value>,
                template<typename> typename PVector = types::vector1d_t, 
                typename KSubVector = typename KVector::value_type,
                typename... Parameters>
            auto solve_p(const AVector& arguments, const VVector& initial_value, const Function& function,
                         const PVector<std::tuple<Parameters...>>& parameters) const 
            {
                return solver<RVector>::solve(
                    create_lazy_bound_p<Function, KVector, AVector, VVector, PVector, KSubVector, Parameters...>(
                        arguments, initial_value, function, parameters
                    )
                );
            }

            template<
                typename Function,
                typename RVector = types::vector2d_t<Value>,
                typename KVector = types::vector2d_t<Value>,
                typename VVector = types::vector1d_t<Argument>, 
                template<typename> typename PVector = types::vector1d_t, 
                typename KSubVector = typename KVector::value_type,
                typename... Parameters>
            auto solve_uniform_p(const Argument& a, const Argument& b, const size_t n, const VVector& initial_value,
                                 const Function& function, const PVector<std::tuple<Parameters...>>& parameters) const 
            {
                return solver<RVector>::solve(
                    create_lazy_uniform_bound_p<Function, KVector, VVector, PVector, KSubVector, Parameters...>(
                        a, b, n, initial_value, function, parameters
                    )
                );
            }

        private:

            const Coeffs _coefficients;


            template<typename RVector>
            struct solver {

                template<typename RK>
                static auto solve(RK&& rk) {
                    auto result = utils::constructor<RVector>::construct(
                        utils::arguments(rk.size()),
                        utils::no_arguments()
                    );
                    assign(result[0], rk.initial_value());
                    for (size_t i = 1; i < rk.size(); ++i)
                        assign(result[i], rk());
                    return result;
                }

            };

            template<typename T1, typename T2>
            static void assign(T1& to, T2&& from) {
                utils::assign<T1, T2, utils::move_assigner, utils::copy_assigner, utils::loop_assigner>(to, std::move(from));
            }

    };

}// namespace ssh
