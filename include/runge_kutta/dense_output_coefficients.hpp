#pragma once

#include "../utils/types.hpp"
#include "../utils/utils.hpp"
#include <functional>

namespace dork {

    template<typename Vector2, typename T, typename V = T>
    class dense_output_coefficients {

    public:

        using functions_t = types::vector1d_t<std::function<V(const T&)>>;
        using argument_t = T;
        using value_t = V;

        dense_output_coefficients() = default;
        ~dense_output_coefficients() = default;

        explicit dense_output_coefficients(const Vector2 coefficients) :
                _coefficients(coefficients) {}

        template<size_t N = 0>
        auto generate_functions() const {
            functions_t result;
            types::vector1d_t<T> coeffs(width(), T(0));
            for (size_t i = N == 0 ? 1 : N; i <= width(); ++i) {
                coeffs[i - 1] = utils::factorial<T>(i) / utils::factorial<T>(i - N);
            }
            auto coefficients = _coefficients;
            for (auto& it : coefficients)
                for (size_t i = 0; i < it.size(); ++i)
                    it[i] *= coeffs[i];
            for (const auto& it : coefficients) {
                result.push_back([it](const T& value) {
                    auto result = V(0);
                    auto v = N == 0 ? value : T(1);
                    for (size_t i = N == 0 ? 0 : N - 1; i < it.size(); ++i) {
                        result += it[i] * v;
                        v *= value;
                    }
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

    template<typename T, typename V = T>
    class runge_kutta_dense_output_coefficients : public dense_output_coefficients<std::array<std::array<T, 3>, 4>, T, V> {

    public:

        static constexpr std::array<std::array<T, 3>, 4> coefficients =
        {
            T(1), -T(3) / T(2),  T(2) / T(3),
            T(0),  T(1),        -T(2) / T(3),
            T(0),  T(1),        -T(2) / T(3),
            T(0), -T(1) / T(2),  T(2) / T(3)
        };

        runge_kutta_dense_output_coefficients() : dense_output_coefficients<std::array<std::array<T, 3>, 4>, T, V>(coefficients) {}

    };

    template<typename T, typename V = T>
    class dormand_prince_dense_output_coefficients1 : public dense_output_coefficients<std::array<std::array<T, 5>, 7>, T, V> {

    public:

        static constexpr std::array<std::array<T, 5>, 7> coefficients =
        {
            T(1), -T(4034104133) / T(1410260304),     T(105330401) / T(33982176),  -T(13107642775) / T(11282082432),    T(6542295) / T(470086768),
            T(0),  T(0),                              T(0),                         T(0),                               T(0),
            T(0),  T(132343189600) / T(32700410799), -T(833316000) / T(131326951),  T(91412856700) / T(32700410799),   -T(523383600) / T(10900136933),
            T(0), -T(115792950) / T(29380423),        T(185270875) / T(16991088),  -T(12653452475) / T(1880347072),     T(98134425) / T(235043384),
            T(0),  T(70805911779) / T(24914598704),  -T(4531260609) / T(600351776), T(988140236175) / T(199316789632), -T(14307999165) / T(24914598704),
            T(0), -T(331320693) / T(205662961),       T(31361737) / T(7433601),    -T(2426908385) / T(822651844),       T(97305120) / T(205662961),
            T(0),  T(44764047) / T(29380423),        -T(1532549) / T(353981),       T(90730570) / T(29380423),         -T(8293050) / T(29380423)
        };

        dormand_prince_dense_output_coefficients1() : dense_output_coefficients<std::array<std::array<T, 5>, 7>, T, V>(coefficients) {}

    };

    template<typename T, typename V = T>
    class dormand_prince_dense_output_coefficients2 : public dense_output_coefficients<std::array<std::array<T, 5>, 7>, T, V> {

    public:

        static constexpr std::array<std::array<T, 5>, 7> coefficients =
        {
                T(1), -T(1615727663233) / T(564104121600),    T(31624252951) / T(10194652800),  -T(13107642775) / T(11282082432),    T(6542295) / T(470086768),
                T(0),  T(0),                                  T(0),                              T(0),                               T(0),
                T(0),  T(663801958033) / T(163502053995),    -T(5364212186) / T(844244685),      T(91412856700) / T(32700410799),   -T(523383600) / T(10900136933),
                T(0), -T(76193498033) / T(18803470720),       T(1243516717) / T(113273920),     -T(12653452475) / T(1880347072),     T(98134425) / T(235043384),
                T(0),  T(29843066025657) / T(9965839481600), -T(459233295093) / T(60035177600),  T(988140236175) / T(199316789632), -T(14307999165) / T(24914598704),
                T(0), -T(8929386631) / T(5141574025),         T(114231227) / T(26548575),       -T(2426908385) / T(822651844),       T(97305120) / T(205662961),
                T(0),  T(44764047) / T(29380423),            -T(1532549) / T(353981),            T(90730570) / T(29380423),         -T(8293050) / T(29380423)
        };

        dormand_prince_dense_output_coefficients2() : dense_output_coefficients<std::array<std::array<T, 5>, 7>, T, V>(coefficients) {}

    };

}// namespace dork
