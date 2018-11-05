#pragma once

#include <functional>
#include <cmath>
#include <type_traits>
#include <algorithm>
#include "../utils/types.hpp"
#include "../utils/constructor.hpp"

namespace ssh {

    namespace implementation {

        template<
            typename Argument, 
            size_t N = 0, 
            typename Value = Argument, 
            typename Functions = types::vector1d_t<std::function<Value(const Argument&)>>>
        class dense_output_implementation {
            
            public:

                dense_output_implementation() = delete;
                ~dense_output_implementation() = default;

                dense_output_implementation(const dense_output_implementation&) = default;
                dense_output_implementation(dense_output_implementation&&) = default;

                dense_output_implementation(const Functions& functions) : _functions(functions) {}

                template<typename Vector = types::vector1d_t<Value>>
                auto interpolate(const Argument& h, const Vector& k) {
                    return interpolate<Vector>(h, k, k[0]);
                }

                template<typename Vector = types::vector1d_t<Value>>
                auto interpolate(const Argument& h, const Vector& k, const Value& y) const {
                    return [functions=_functions, h, k, y](const Argument& x) {
                        auto result = Value(0);
                        for (size_t i = 0; i < k.size(); ++i)
                            result += functions[i](x) * k[i];
                        return (N == 0 ? y : Value(0)) + h * result;
                    };
                }

                template<typename KVector = types::vector2d_t<Value>, typename RVector = typename KVector::value_type>
                auto interpolate_p(const Argument& h, const KVector& k) const {
                    return interpolate_p<KVector, RVector>(h, k, k[0]);
                }

                template<typename KVector = types::vector2d_t<Value>, typename RVector = types::vector1d_t<Value>>
                auto interpolate_p(const Argument& h, const KVector& k, const RVector& y) const {
                    return [functions=_functions, h, k, y](const Argument& x) {
                        auto result = utils::constructor<RVector>::construct(
                            utils::arguments(y.size(), Value(0)),
                            utils::arguments(y.size()),
                            utils::no_arguments()
                        );
                        for (size_t i = 0; i < y.size(); ++i) {
                            for (size_t j = 0; j < functions.size(); ++j)
                                result[i] += functions[j](x) * k[i][j];
                            result[i] = (N == 0 ? y[i] : Value(0)) + h * result[i];
                        }
                        return result;
                    };
                }

                template<typename KVector = types::vector2d_t<Value>>
                auto interpolate_bound(const Argument& h, const KVector& k) {
                    return interpolate_bound<KVector>(h, k, k[0]);
                }

                template<typename KVector = types::vector2d_t<Value>>
                auto interpolate_bound(const Argument& h, const KVector& k, const Value& y) {
                    return [functions=_functions, &h, &k, &y](const Argument& x) {
                        auto result = Value(0);
                        for (size_t i = 0; i < k.size(); ++i)
                            result += functions[i](x) * k[i];
                        return (N == 0 ? y : Value(0)) + h * result;
                    };
                }

                template<typename KVector = types::vector2d_t<Value>, typename RVector = typename KVector::value_type>
                auto interpolate_bound_p(const Argument& h, const KVector& k) {
                    return interpolate_bound_p<KVector, RVector>(h, k, k[0]);
                }

                template<typename KVector = types::vector2d_t<Value>, typename RVector = types::vector1d_t<Value>>
                auto interpolate_bound_p(const Argument& h, const KVector& k, const RVector& y) {
                    return [functions=_functions, &h, &k, &y](const Argument& x) {
                        auto result = utils::constructor<RVector>::construct(
                            utils::arguments(y.size(), Value(0)),
                            utils::arguments(y.size()),
                            utils::no_arguments()
                        );
                        for (size_t i = 0; i < k.size(); ++i) {
                            for (size_t j = 0; j < functions.size(); ++j)
                                result[i] += functions[j](x) * k[i][j];
                            result[i] = (N == 0 ? y[i] : Value(0)) + h * result[i];
                        }
                        return result;
                    };
                }

                template<typename RK>
                auto interpolate_solution(RK& rk) const {
                    return _interpolate_solution<Value, function>(rk, interpolate_lambda, rk.arguments());
                }

                template<typename RK>
                auto interpolate_solution_uniform(RK& rk) const {
                    return _interpolate_solution<Value, uniform_function>(rk, interpolate_lambda, rk.bounds().first, rk.step());
                }

                template<typename RK, typename VVector = typename RK::vvector_t>
                auto interpolate_solution_p(RK& rk) const {
                    return _interpolate_solution<VVector, function>(rk, interpolate_p_lambda, rk.arguments());
                }

                template<typename RK, typename VVector = typename RK::vvector_t>
                auto interpolate_solution_uniform_p(RK& rk) const {
                    return _interpolate_solution<VVector, uniform_function>(rk, interpolate_p_lambda, rk.bounds().first, rk.step());
                }

                template<typename RK, typename AVector, typename RVector = types::vector1d_t<Value>>
                auto solve(RK& rk, const AVector& arguments) const {
                    return _solve<RVector, init_function>(interpolate_solution(rk), arguments.size(), arguments);
                }

                template<typename RK, typename RVector = types::vector1d_t<Value>>
                auto solve_uniform(RK& rk, const size_t n) const {
                    const auto bounds = rk.bounds();
                    return _solve<RVector, init_function_uniform>(
                        interpolate_solution_uniform(rk), n, bounds.first, (bounds.second - bounds.first) / (n - 1));
                }

                template<typename RK, typename AVector, typename RVector = types::vector1d_t<typename RK::vvector_t>>
                auto solve_p(RK& rk, const AVector& arguments) const {
                    return _solve<RVector, init_function>(interpolate_solution_p(rk), arguments.size(), arguments);
                }

                template<typename RK, typename RVector = types::vector1d_t<typename RK::vvector_t>>
                auto solve_uniform_p(RK& rk, const size_t n) const {
                    const auto bounds = rk.bounds();
                    return _solve<RVector, init_function_uniform>(
                        interpolate_solution_uniform_p(rk), n, bounds.first, (bounds.second - bounds.first) / (n - 1));
                }

            private:

                const Functions _functions;

                static constexpr auto interpolate_lambda = 
                    [](const dense_output_implementation* obj, const Argument& h, const auto& k, const Value& y) { 
                            return obj->interpolate(h, k, y);
                    };

                static constexpr auto interpolate_p_lambda = 
                    [](const dense_output_implementation* obj, const Argument& h, const auto& k, const auto& y) { 
                            return obj->interpolate_p(h, k, y);
                    };

                template<typename Fs, typename As>
                class function {

                    public:

                        function() = delete;
                        ~function() = default;

                        function(const function&) = default;
                        function(function&&) = default;

                        function(Fs&& functions, const As& arguments) : 
                            _functions(std::move(functions)), _arguments(arguments) {}

                        auto operator()(const Argument& x) const {
                            size_t index;
                            if (x <= _arguments.front())
                                index = 0;
                            else if (x >= _arguments.back())
                                index = _arguments.size() - 2;
                            else
                                index = std::distance(_arguments.begin(), 
                                    std::upper_bound(_arguments.begin(), _arguments.end(), x)) - 1;
                            const auto a = _arguments[index];
                            return _functions[index]((x - a) / (_arguments[index + 1] - a));
                        }

                        auto domain() const {
                            return std::make_pair(_arguments.front(), _arguments.back());
                        }
                    
                    private:

                        const Fs _functions;
                        const As _arguments;

                };

                template<typename Fs, typename = void>
                class uniform_function {

                    public:

                        uniform_function() = delete;
                        ~uniform_function() = default;

                        uniform_function(const uniform_function&) = default;
                        uniform_function(uniform_function&&) = default;

                        uniform_function(Fs&& functions, const Argument& a, const Argument& d) : 
                            _functions(std::move(functions)), _a(a), _d(d) {}

                        auto operator()(const Argument& x) const {
                            const auto index = std::min<size_t>(
                                static_cast<size_t>(std::max(Argument(0), std::floor((x - _a) / _d))), 
                                _functions.size() - 1
                            );
                            return _functions[index]((x - _a - _d * index) / _d);
                        }

                        auto domain() const {
                            return std::make_pair(_a, _a + _functions.size() * _d);
                        }
                    
                    private:

                        const Fs _functions;
                        const Argument _a, _d;

                };

                template<typename F, typename A>
                class init_function {

                public:

                    init_function(const F& f, const A& arguments) : _f(f), _arguments(arguments) {}

                    auto operator()(const size_t i) const {
                        return _f(_arguments[i]);
                    }

                private:

                    const F _f;
                    const A& _arguments;

                };

                template<typename F, typename = void>
                class init_function_uniform {

                public:

                    init_function_uniform(const F& f, const Argument& a, const Argument& d) : _f(f), _a(a), _d(d) {}

                    auto operator()(const size_t i) const {
                        return _f(_a + i * _d);
                    }

                private:

                    const F _f;
                    const Argument _a, _d;

                };

                template<typename RK, typename V, typename Interpolate>
                auto generate_interpolating_functions(RK& rk, const Interpolate& function) const {
                    auto k = rk.construct_kvector();
                    types::vector1d_t<std::function<V(const Argument&)>> functions(rk.size());
                    auto y0 = rk.initial_value();
                    for (size_t i = 0; i < rk.size(); ++i) {
                        auto y = rk(k);
                        functions[i] = function(this, rk.step(i), k, y0);
                        y0 = y;
                    }
                    return functions;
                }

                template<typename V, template<typename, typename> typename F, typename RK, typename IF, typename... Args>
                auto _interpolate_solution(RK& rk, const IF& function, const Args&... args) const {
                    return F(generate_interpolating_functions<RK, V>(rk, function), args...);
                }

                template<typename RVector, template<typename, typename> typename IF, typename F, typename... Args>
                auto _solve(const F& f, const size_t n, const Args&... args) const {
                    auto result = utils::constructor<RVector>::construct(
                        utils::arguments(n),
                        utils::no_arguments()
                    );
                    utils::init_vector(result, IF(f, args...));
                    return result;
                }

        };

    }// namespace implementation

    template<size_t N>
    struct dense_output_derivative {

        template<
            typename Coeffs,
            typename Argument = typename Coeffs::argument_t,
            typename Value = typename Coeffs::value_t>
        static auto from_coefficients() { 
            return implementation::dense_output_implementation<Argument, N, Value>(Coeffs().template generate_functions<N>());
        }

        template<
            typename Coeffs,
            typename Argument = typename Coeffs::argument_t,
            typename Value = typename Coeffs::value_t>
        static auto from_coefficients(const Coeffs& coefficients) {
            return implementation::dense_output_implementation<Argument, N, Value>(coefficients.template generate_functions<N>());
        }

        template<
            template<typename> typename Coeffs,
            typename Argument,
            typename Value = typename Coeffs<Argument>::value_t>
        static auto from_coefficients(const Coeffs<Argument>& coefficients) {
            return implementation::dense_output_implementation<Argument, N, Value>(coefficients.template generate_functions<N>());
        }

        template<
            typename Argument,
            typename Value = Argument,
            typename Functions = types::vector1d_t<std::function<Value(const Argument&)>>>
        static auto from_functions(const Functions& functions) {
            return implementation::dense_output_implementation<Argument, N, Value>(functions);
        }

    };

    using dense_output = dense_output_derivative<0>;

}// namespace ssh
