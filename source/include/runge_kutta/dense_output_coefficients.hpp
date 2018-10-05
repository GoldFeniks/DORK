#pragma once

#include "../utils/types.hpp"
#include "../utils/utils.hpp"
#include <functional>

namespace ssh {

    template<typename Vector2, typename T>
    class dense_output_coefficients {

    public:

        dense_output_coefficients() = default;
        ~dense_output_coefficients() = default;

        dense_output_coefficients(const Vector2 coefficients) :
                _coefficients(coefficients) {}

        template<size_t N = 0>
        auto generate_functions() const {
            types::vector1d_t<std::function<T(const T&)>> result;
            types::vector1d_t<T> coeffs(width(), T(0));
            for (size_t i = N == 0 ? 1 : N; i <= coeffs.size(); ++i) {
                coeffs[i - 1] = utils::factorial<T>(i) / utils::factorial<T>(i - N);
            }
            auto coefficients = _coefficients;
            for (auto& it : coefficients)
                for (size_t i = 0; i < it.size(); ++i)
                    it[i] *= coeffs[i];
            for (const auto& it : coefficients) {
                result.push_back([it](const T& value) {
                    auto result = T(0);
                    auto v = T(1);
                    for (const auto& c : it)
                        result += c * (v *= value);
                    return result;
                });
            }
            return result;
        }

        const T& get(const size_t i, const size_t j) {
            return _coefficients[i][j];
        }

        const size_t height() const {
            return _coefficients.size();
        }

        const size_t width() const {
            return _coefficients[0].size();
        }

    private:

        const Vector2 _coefficients;

    };

    template<typename T, size_t N>
    class runge_kutta_dense_output_coefficients {};

    template<typename T>
    class runge_kutta_dense_output_coefficients<T, 3> : public dense_output_coefficients<std::array<std::array<T, 3>, 4>, T> {

    public:

        static constexpr std::array<std::array<T, 3>, 4> coefficients =
            {
                T(1), -T(3) / T(2),  T(2) / T(3),
                T(0),  T(1),        -T(2) / T(3),
                T(0),  T(1),        -T(2) / T(3),
                T(0), -T(1) / T(2),  T(2) / T(3)
            };

        runge_kutta_dense_output_coefficients() : dense_output_coefficients<std::array<std::array<T, 3>, 4>, T>(coefficients) {}

    };

}// namespace ssh
