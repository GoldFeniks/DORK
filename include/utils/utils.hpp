#pragma once

#include <tuple>
#include <vector>
#include <cstddef>
#include <functional>
#include <type_traits>
#include "types.hpp"

namespace dork {

    namespace utils {

        template<typename...>
        struct collect_types { collect_types() = delete; };

        template<typename T, typename... V>
        struct collect_types<T, V...> {

            collect_types() = delete;

            using type = typename collect_types<T, typename collect_types<V...>::type>::type;

        };

        template<typename... T, typename... V>
        struct collect_types<std::tuple<T...>, V...> {

            collect_types() = delete;

            using type = typename collect_types<std::tuple<T...>, typename collect_types<V...>::type>::type;

        };

        template<typename... T, typename... V>
        struct collect_types<std::tuple<T...>, std::tuple<V...>> {

            collect_types() = delete;

            using type = std::tuple<T..., V...>;

        };

        template<typename T, typename... V>
        struct collect_types<T, std::tuple<V...>> {

            collect_types() = delete;

            using type = std::tuple<T, V...>;

        };

        template<typename... T>
        struct collect_types<std::tuple<T...>> {

            collect_types() = delete;

            using type = std::tuple<T...>;

        };

        template<typename T>
        struct collect_types<T> {

            collect_types() = delete;

            using type = std::tuple<T>;

        };

        template<typename... T>
        using collect_types_t = typename collect_types<T...>::type;

        template<typename... T>
        struct const_reference_tuple;

        template<typename T, typename... V>
        struct const_reference_tuple<T, V...> {

            const_reference_tuple() = delete;

            using type = collect_types_t<const T&, typename const_reference_tuple<V...>::type>;

            static type construct(const T& value, const V&... values) {
                return std::tuple_cat(std::tuple<const T&>(value), const_reference_tuple<V...>::construct(values...));
            }

        };

        template<typename... T, typename... V>
        struct const_reference_tuple<std::tuple<T...>, V...> {

            const_reference_tuple() = delete;

            using type = collect_types_t<const T&..., typename const_reference_tuple<V...>::type>;

            static type construct(const std::tuple<T...>& value, const V&... values) {
                return std::tuple_cat(std::apply(
                        [](auto&... values) {
                            return std::tuple<const T&...>(values...);
                        },
                        value
                    ),
                    const_reference_tuple<V...>::construct(values...)
                );
            }

        };

        template<typename... T>
        struct const_reference_tuple<std::tuple<T...>> {

            const_reference_tuple() = delete;

            using type = std::tuple<const T&...>;

            static type construct(const std::tuple<T...>& value) {
                return std::apply([](const auto&... values) { return std::tuple<const T&...>(values...); }, value);
            }

        };

        template<typename T>
        struct const_reference_tuple<T> {

            const_reference_tuple() = delete;

            using type = std::tuple<const T&>;

            static type construct(const T& value) {
                return std::tuple<const T&>(value);
            }

        };

        template<typename... T>
        using const_reference_tuple_t = typename const_reference_tuple<T...>::type;

        template<typename... T>
        auto make_const_reference_tuple(const T&... values) {
            return const_reference_tuple<T...>::construct(values...);
        }

    }// namespace utils

}// namespace dork
