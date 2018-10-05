#pragma once

#include <functional>
#include "../utils/types.hpp"

namespace ssh {

    class dense_output {

    public:

        template<typename T, typename Coeffs, typename Vector = types::vector1d_t<T>, size_t N = 0>
        static auto create(const Coeffs& coefficients, const T& h, const T& y0, const Vector& k) {
            const auto functions = coefficients.template generate_functions<N>();
            return [functions, h, y0, k](const T x) {
                auto result = T(0);
                for (size_t i = 0; i < k.size(); ++i)
                    result += functions[i](x) * k[i];
                return (N == 0 ? y0 : 0) + h * result;
            };
        }

        template<typename T, typename Coeffs, typename Vector = types::vector1d_t<T>, size_t N = 0>
        static auto create_bound(const Coeffs& coefficients, const T& h, const T& y0, const Vector& k) {
            const auto functions = coefficients.template generate_functions<N>();
            return [functions, h, &y0, &k](const T x) {
                auto result = T(0);
                for (size_t i = 0; i < k.size(); ++i)
                    result += functions[i](x) * k[i];
                return (N == 0 ? y0 : 0) + h * result;
            };
        }

    };

    class runge_kutta_dense_output {

    public:

        template<typename T, typename RK, typename Coeffs, typename Vector = types::vector1d_t<T>>
        static Vector solve(RK& rk, const T& h, const T& y0, const size_t n, const Coeffs& coefficients) {
            types::vector1d_t<T> k(coefficients.height());
            auto v = y0;
            Vector result;
            auto output = dense_output::create_bound(coefficients, h, v, k);
            while (rk) {
                auto v1 = rk(k);
                for (size_t i = 0; i < n; ++i)
                    result.push_back(output(T(1) / T(n) * i));
                v = v1;
            }
            result.push_back(v);
            return result;
        }

    };



}// namespace ssh
