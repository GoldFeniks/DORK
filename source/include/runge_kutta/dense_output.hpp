#pragma once

#include <functional>
#include <cmath>
#include <type_traits>
#include "../utils/types.hpp"

namespace ssh {

    template<typename Argument, size_t N = 0, typename Value = Argument, typename Functions = types::vector1d_t<std::function<Value(const Argument&)>>>
    class dense_output {
        
        public:

            dense_output() = delete;
            ~dense_output() = default;

            dense_output(const dense_output&) = default;
            dense_output(dense_output&&) = default;

            dense_output(const Functions& functions) : _functions(functions) {}

            template<typename Vector = types::vector1d_t<Value>>
            auto interpolate(const Argument& h, const Vector& k) {
                return interpolate<Vector>(h, k, k[0]);
            }

            template<typename Vector = types::vector1d_t<Value>>
            auto interpolate(const Argument& h, const Vector& k, const Value& y) {
                return [functions=_functions, h, k, y](const Argument& x) {
                    auto result = Value(0);
                    for (size_t i = 0; i < k.size(); ++i)
                        result += functions[i](x) * k[i];
                    return (N == 0 ? y : Value(0)) + h * result * std::pow(2, N);
                };
            }

            template<typename KVector = types::vector2d_t<Value>, typename RVector = typename KVector::value_type>
            auto interpolate_p(const Argument& h, const KVector& k) {
                return interpolate_p<KVector, RVector>(h, k, k[0]);
            }

            template<typename KVector = types::vector2d_t<Value>, typename RVector = types::vector1d_t<Value>>
            auto interpolate_p(const Argument& h, const KVector& k, const RVector& y) {
                return [functions=_functions, h, k, y](const Argument& x) {
                    auto result = std::is_constructible_v<RVector, size_t, Value>() ? RVector(k.size(), Value(0)) : RVector();
                    for (size_t i = 0; i < k.size(); ++i) {
                        for (size_t j = 0; j < functions.size(); ++j)
                            result[i] += functions[j](x) * k[i][j];
                        result[i] = (N == 0 ? y[i] : Value(0)) + h * result[i] * std::pow(2, N);
                    }
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
                    return (N == 0 ? y : Value(0)) + h * result * std::pow(2, N);
                };
            }

            template<typename KVector = types::vector2d_t<Value>, typename RVector = typename KVector::value_type>
            auto interpolate_bound_p(const Argument& h, const KVector& k) {
                return interpolate_bound_p<KVector, RVector>(h, k, k[0]);
            }

            template<typename KVector = types::vector2d_t<Value>, typename RVector = types::vector1d_t<Value>>
            auto interpolate_bound_p(const Argument& h, const KVector& k, const RVector& y) {
                return [functions=_functions, &h, &k, &y](const Argument& x) {
                    auto result = std::is_constructible_v<RVector, size_t, Value>() ? RVector(k.size(), Value(0)) : RVector();
                    for (size_t i = 0; i < k.size(); ++i) {
                        for (size_t j = 0; j < functions.size(); ++j)
                            result[i] += functions[j](x) * k[i][j];
                        result[i] = (N == 0 ? y[i] : Value(0)) + h * result[i] * std::pow(2, N);
                    }
                };
            }

            template<typename Coeffs>
            static auto from_coefficients() { 
                return dense_output(Coeffs().template generate_functions<N>());
            }

            template<typename Coeffs>
            static auto from_coefficients(const Coeffs& coefficients) {
                return dense_output(coefficients.template generate_functions<N>());
            }

            template<typename RK, typename RVector = types::vector1d_t<Value>>
            auto solve(RK& rk, const Argument& h, const Value& y0, const size_t n) {
                types::vector1d_t<Value> k(_functions.size());
                auto v = y0;
                RVector result;
                auto output = interpolate_bound(h, k, v);
                while (rk) {
                    result.push_back(v);
                    auto v1 = rk(k);
                    for (size_t i = 1; i < n; ++i)
                        result.push_back(output(Argument(1) / Argument(n) * Argument(i)));
                    v = v1;
                }
                result.push_back(v);
                return result;
            }

        private:

            const Functions _functions;

    };

}// namespace ssh
